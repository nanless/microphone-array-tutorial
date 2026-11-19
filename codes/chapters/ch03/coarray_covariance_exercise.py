"""E03-07: reconstruct small co-array statistics from physical covariances.

Run with ``.venv/bin/python -m codes.chapters.ch03.coarray_covariance_exercise``.
The numbers are exact narrowband teaching inputs, not recorded microphone data.
The finite-snapshot example shows why a lag-filled Toeplitz matrix must be
checked before it is called a usable covariance matrix.
"""

from __future__ import annotations

import json

import numpy as np


POSITIONS = np.array([0, 1, 3], dtype=int)  # half-wavelength grid units


def _is_hermitian(matrix: np.ndarray) -> bool:
    """Compare at a common safe scale; zero is exactly Hermitian."""
    scale = max(float(np.max(abs(matrix.real))), float(np.max(abs(matrix.imag))))
    if scale == 0:
        return True
    scaled = matrix.real / scale + 1j * (matrix.imag / scale)
    return bool(np.allclose(scaled, scaled.conj().T, atol=1e-12, rtol=1e-12))


def _contains_bool(value: object) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return True
    if isinstance(value, np.ndarray):
        return value.dtype.kind == "b"
    if isinstance(value, (list, tuple)):
        return any(_contains_bool(item) for item in value)
    return False


def average_ordered_lags(covariance: np.ndarray, positions: np.ndarray = POSITIONS) -> dict[int, complex]:
    """Average finite Hermitian R[i,j] over signed integer p[i]-p[j].

    This is a reordering tool, not a PSD validator. The caller must establish
    that an input described as a physical covariance has the required model.
    Python integers preserve differences exceeding a fixed-width integer dtype.
    """
    raw_cov = np.asarray(covariance)
    coordinates = np.asarray(positions)
    if (raw_cov.dtype.kind not in "iufc" or coordinates.dtype.kind not in "iuf"
            or _contains_bool(covariance) or _contains_bool(positions)):
        raise ValueError("require numeric non-bool covariance and integer positions")
    matrix = np.asarray(raw_cov, dtype=complex)
    if (matrix.shape != (coordinates.size, coordinates.size)
            or coordinates.ndim != 1 or coordinates.size == 0
            or np.any(~np.isfinite(matrix)) or np.any(~np.isfinite(coordinates))
            or (coordinates.dtype.kind == "f" and np.any(coordinates != np.rint(coordinates)))):
        raise ValueError("require finite square covariance and integer one-dimensional positions")
    if not _is_hermitian(matrix):
        raise ValueError("input statistics must be Hermitian")
    coordinates = [int(p) for p in coordinates]
    grouped: dict[int, list[complex]] = {}
    for i, pi in enumerate(coordinates):
        for j, pj in enumerate(coordinates):
            grouped.setdefault(pi-pj, []).append(matrix[i, j])
    result = {}
    for lag, values in sorted(grouped.items()):
        # Preserve ordinary finite means. Only an overflowing accumulation
        # needs division before summing; a finite average can still exist.
        with np.errstate(over="ignore", invalid="ignore"):
            mean = np.mean(values)
        if not np.isfinite(mean):
            mean = sum(value / len(values) for value in values)
        if not np.isfinite(mean):
            raise ValueError("lag average exceeds supported floating-point range")
        result[lag] = complex(mean)
    return result


def virtual_toeplitz(lags: dict[int, complex], size: int) -> np.ndarray:
    """Fill a contiguous virtual ULA; this operation alone does not ensure PSD."""
    if isinstance(size, (bool, np.bool_)) or not isinstance(size, (int, np.integer)) or int(size) < 1:
        raise ValueError("size must be a positive integer")
    size_int = int(size)
    if not isinstance(lags, dict):
        raise ValueError("lags must be a dict mapping integer lags to finite numeric scalars")
    normalized_lags: dict[int, complex] = {}
    for raw_lag, raw_val in lags.items():
        if isinstance(raw_lag, (bool, np.bool_)) or not isinstance(raw_lag, (int, np.integer)):
            raise ValueError("lag keys must be integers, not bool")
        val_arr = np.asarray(raw_val)
        if isinstance(raw_val, (bool, np.bool_)) or val_arr.ndim != 0 or val_arr.dtype.kind not in "iufc":
            raise ValueError("lag values must be numeric scalars, not bool")
        normalized_lags[int(raw_lag)] = complex(val_arr)
    required = set(range(1 - size_int, size_int))
    if not required.issubset(normalized_lags):
        raise ValueError("all signed lags in the requested contiguous segment are required")
    result = np.array(
        [[normalized_lags[i - j] for j in range(size_int)] for i in range(size_int)],
        dtype=complex,
    )
    if not np.all(np.isfinite(result)):
        raise ValueError("lag values must be finite")
    if not _is_hermitian(result):
        raise ValueError("lag values must satisfy conjugate symmetry")
    return result


def run_exercise() -> dict:
    """Return an ideal PSD case and a finite-snapshot indefinite counterexample."""
    # Broadside z=1 and +30-degree z=i: a_m = z**p_m under this book's
    # positive-angle, half-wavelength ULA convention. Each has unit power.
    a_broadside = np.ones(3, dtype=complex)
    a_30 = (1j) ** POSITIONS
    physical = (np.outer(a_broadside, a_broadside.conj())
                + np.outer(a_30, a_30.conj()) + 0.1*np.eye(3))
    lags = average_ordered_lags(physical)
    virtual = virtual_toeplitz(lags, 4)

    # A single finite sample is a valid rank-one physical covariance.
    # Its lag-filled virtual Toeplitz matrix has a negative eigenvalue.
    sample = np.array([1, 0, 1], dtype=complex)
    sample_physical = np.outer(sample, sample.conj())
    sample_lags = average_ordered_lags(sample_physical)
    sample_virtual = virtual_toeplitz(sample_lags, 4)

    encode = lambda values: {str(k): [float(v.real), float(v.imag)] for k, v in values.items()}
    return {
        "physical_positions_half_wavelength": POSITIONS.tolist(),
        "physical_channels": 3,
        "virtual_positions_half_wavelength": [0, 1, 2, 3],
        "source_directions_deg": [0, 30],
        "source_powers": [1.0, 1.0],
        "independent_channel_noise_power": 0.1,
        "ideal_lags_real_imag": encode(lags),
        "ideal_virtual_eigenvalues": np.linalg.eigvalsh(virtual).tolist(),
        "sample_snapshot_real_imag": [[float(z.real), float(z.imag)] for z in sample],
        "sample_lags_real_imag": encode(sample_lags),
        "sample_virtual_eigenvalues": np.linalg.eigvalsh(sample_virtual).tolist(),
        "sample_physical_eigenvalues": np.linalg.eigvalsh(sample_physical).tolist(),
    }


if __name__ == "__main__":
    print(json.dumps(run_exercise(), ensure_ascii=False, indent=2))
