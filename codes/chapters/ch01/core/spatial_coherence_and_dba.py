"""E01-11 & E01-12: spatial noise coloring, IEC 61672-1 A-weighting, and SOFA coordinates.

Pure-NumPy deterministic kernels for Chapter 1 without Chapter 2+ imports or SciPy.
No files, directories, experiments, or network operations run at import time.
"""

from __future__ import annotations

import math

import numpy as np


IEC61672_POLE_FREQUENCIES_HZ: tuple[float, float, float, float] = (
    20.598997,
    107.65265,
    737.86223,
    12194.217,
)


def _finite_real_scalar(value: object, name: str) -> float:
    """Validate one finite real scalar, explicitly rejecting booleans and arrays."""
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a real number, not bool")
    raw = np.asarray(value)
    if raw.ndim != 0 or raw.dtype.kind not in "iuf":
        raise ValueError(f"{name} must be a finite real scalar")
    with np.errstate(over="ignore", invalid="ignore"):
        scalar = float(raw)
    if not math.isfinite(scalar):
        raise ValueError(f"{name} must be finite and representable as float64")
    return scalar


def _finite_real_array(value: object, name: str) -> np.ndarray:
    """Validate finite real array input without discarding imaginary parts."""
    raw = np.asarray(value)
    if raw.dtype.kind not in "iuf":
        raise ValueError(f"{name} must contain real numeric values, not complex, bool, or text")
    with np.errstate(over="ignore", invalid="ignore"):
        arr = np.asarray(raw, dtype=np.float64)
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be finite and representable as float64")
    return arr


def hermitian_psd_sqrt_coloring(
    coherence: np.ndarray,
    *,
    hermitian_rtol: float = 1e-10,
    psd_rtol: float = 1e-10,
) -> dict[str, object]:
    """Compute unique Hermitian PSD square-root coloring matrix U @ Lambda^(1/2) @ U^H.

    Validates finite 2-D square numeric input, checks relative Hermitian symmetry
    within ``hermitian_rtol``, computes ``np.linalg.eigh`` on the symmetrized matrix,
    rejects genuinely negative eigenvalues ``< -psd_rtol * max_scale`` with
    ``np.linalg.LinAlgError``, clamps near-zero eigenvalues in
    ``[-psd_rtol * max_scale, psd_rtol * max_scale]`` to ``0.0``, and returns eigenvalues, numerical
    rank, Hermitian square-root matrix, and reconstruction error.
    """
    h_tol = _finite_real_scalar(hermitian_rtol, "hermitian_rtol")
    p_tol = _finite_real_scalar(psd_rtol, "psd_rtol")
    if h_tol < 0.0 or p_tol < 0.0:
        raise ValueError("hermitian_rtol and psd_rtol must be non-negative")

    raw = np.asarray(coherence)
    if raw.dtype.kind not in "iufc":
        raise ValueError("coherence must be a numeric 2-D square matrix, not bool or text")
    is_real_input = raw.dtype.kind in "iuf"
    with np.errstate(over="ignore", invalid="ignore"):
        mat = np.asarray(raw, dtype=np.float64 if is_real_input else np.complex128)
    if mat.ndim != 2 or mat.shape[0] < 1 or mat.shape[0] != mat.shape[1] or not np.all(np.isfinite(mat)):
        raise ValueError("coherence must be a non-empty finite 2-D square matrix")

    peak = float(np.max(np.abs(mat) if is_real_input else np.maximum(np.abs(mat.real), np.abs(mat.imag))))
    if peak == 0.0:
        zeros = np.zeros_like(mat)
        return {
            "num_channels": int(mat.shape[0]),
            "is_real": bool(is_real_input),
            "eigenvalues": np.zeros(mat.shape[0], dtype=np.float64),
            "numerical_rank": 0,
            "coloring_matrix": zeros,
            "reconstruction_max_abs_error": 0.0,
        }

    reduced = mat / peak if is_real_input else (mat.real / peak + 1j * (mat.imag / peak))
    asym = float(np.max(np.abs(reduced - reduced.conj().T)))
    if asym > h_tol:
        raise ValueError(
            f"coherence must be Hermitian within relative tolerance {h_tol:g} (got {asym:.3e})"
        )

    sym_reduced = 0.5 * reduced + 0.5 * reduced.conj().T
    sym = 0.5 * mat + 0.5 * mat.conj().T
    evals_red, evecs = np.linalg.eigh(sym_reduced)
    max_scale = max(1.0, float(np.max(np.abs(evals_red))))
    if float(evals_red[0]) < -p_tol * max_scale:
        raise np.linalg.LinAlgError(
            f"coherence matrix is not positive semidefinite (min scaled eigenvalue {evals_red[0]:.6e})"
        )

    evals_clamped_red = np.where(evals_red > p_tol * max_scale, evals_red, 0.0)
    rank = int(np.sum(evals_red > p_tol * max_scale))
    with np.errstate(over="ignore", invalid="ignore"):
        evals = evals_clamped_red * peak
        sqrt_evals = np.sqrt(evals_clamped_red) * math.sqrt(peak)
        sqrt_mat = (evecs * sqrt_evals[np.newaxis, :]) @ evecs.conj().T
        sqrt_mat = 0.5 * sqrt_mat + 0.5 * sqrt_mat.conj().T
        if is_real_input:
            sqrt_mat = np.asarray(sqrt_mat.real, dtype=np.float64)
        recon = sqrt_mat @ sqrt_mat.conj().T
        recon_err = float(np.max(np.abs(recon - sym)))
    if not (np.all(np.isfinite(evals)) and np.all(np.isfinite(sqrt_mat)) and math.isfinite(recon_err)):
        raise ValueError("PSD square-root coloring exceeded float64 range")

    return {
        "num_channels": int(mat.shape[0]),
        "is_real": bool(is_real_input),
        "eigenvalues": evals,
        "numerical_rank": rank,
        "coloring_matrix": sqrt_mat,
        "reconstruction_max_abs_error": recon_err,
    }


def equicorrelated_coherence(
    num_channels: int,
    rho: float,
) -> np.ndarray:
    """Construct M x M unit-diagonal equicorrelated matrix (1 - rho)*I + rho*1*1^T.

    Validates integer ``num_channels >= 2`` (rejecting bool) and finite real scalar ``rho``.
    Does NOT clip ``rho`` so callers can test PSD rejection below ``-1/(M - 1)``.
    """
    if isinstance(num_channels, (bool, np.bool_)) or not isinstance(num_channels, (int, np.integer)) or num_channels < 2:
        raise ValueError("num_channels must be an integer >= 2")
    rho_val = _finite_real_scalar(rho, "rho")
    m = int(num_channels)
    mat = np.full((m, m), rho_val, dtype=np.float64)
    np.fill_diagonal(mat, 1.0)
    return mat


def color_white_samples(
    coloring_matrix: np.ndarray,
    white_samples: np.ndarray,
) -> np.ndarray:
    """Apply M x M coloring matrix C to M x N white samples Z -> V = C @ Z."""
    c_raw = np.asarray(coloring_matrix)
    z_raw = np.asarray(white_samples)
    if c_raw.dtype.kind not in "iufc" or z_raw.dtype.kind not in "iufc":
        raise ValueError("coloring_matrix and white_samples must be numeric, not bool or text")
    is_real = c_raw.dtype.kind in "iuf" and z_raw.dtype.kind in "iuf"
    dtype = np.float64 if is_real else np.complex128
    with np.errstate(over="ignore", invalid="ignore"):
        c_mat = np.asarray(c_raw, dtype=dtype)
        z_mat = np.asarray(z_raw, dtype=dtype)
    if c_mat.ndim != 2 or c_mat.shape[0] < 1 or c_mat.shape[0] != c_mat.shape[1] or not np.all(np.isfinite(c_mat)):
        raise ValueError("coloring_matrix must be a non-empty finite square matrix")
    if z_mat.ndim != 2 or z_mat.shape[0] != c_mat.shape[0] or z_mat.shape[1] < 1 or not np.all(np.isfinite(z_mat)):
        raise ValueError("white_samples must be finite with shape (M, N) matching coloring_matrix")
    with np.errstate(over="ignore", invalid="ignore"):
        colored = c_mat @ z_mat
    if not np.all(np.isfinite(colored)):
        raise ValueError("colored output exceeded float64 range")
    return colored


def _r_a_raw(freq: np.ndarray) -> np.ndarray:
    """Evaluate unnormalized IEC 61672-1 R_A(f) safely without f^4 overflow."""
    f1, f2, f3, f4 = IEC61672_POLE_FREQUENCIES_HZ
    f_sq = freq * freq
    # Factor R_A(f) = [f^2 / (f^2 + f1^2)] * [f4^2 / (f^2 + f4^2)] * [f / sqrt(f^2 + f2^2)] * [f / sqrt(f^2 + f3^2)]
    # Each factor lies strictly in (0, 1), preventing intermediate f^4 overflow for high frequencies.
    term1 = f_sq / (f_sq + f1 * f1)
    term4 = (f4 * f4) / (f_sq + f4 * f4)
    term2 = freq / np.sqrt(f_sq + f2 * f2)
    term3 = freq / np.sqrt(f_sq + f3 * f3)
    return term1 * term4 * term2 * term3


def iec61672_a_weighting(frequencies_hz: np.ndarray | float) -> dict[str, np.ndarray | float]:
    """Evaluate IEC 61672-1 relative linear gain, power weight, and dB transfer A(f).

    Validates finite strictly positive real frequencies (Hz), normalizes exactly at
    1000.0 Hz so A(1000 Hz) == 0.0 dB, and guards intermediate powers against overflow.
    """
    raw = np.asarray(frequencies_hz)
    is_scalar = raw.ndim == 0
    if is_scalar:
        f_val = _finite_real_scalar(frequencies_hz, "frequencies_hz")
        freq = np.asarray([f_val], dtype=np.float64)
    else:
        freq = _finite_real_array(frequencies_hz, "frequencies_hz")
        if freq.ndim != 1 or freq.size < 1:
            raise ValueError("frequencies_hz must be a scalar or non-empty 1-D real array")
    if np.any(freq <= 0.0):
        raise ValueError("frequencies_hz must be strictly positive for IEC 61672-1 dB evaluation")

    with np.errstate(over="ignore", invalid="ignore"):
        ra_f = _r_a_raw(freq)
        ra_1000 = float(_r_a_raw(np.asarray([1000.0], dtype=np.float64))[0])
        gain = ra_f / ra_1000
        power_weight = gain * gain
        a_db = 20.0 * np.log10(gain)
    if not (np.all(np.isfinite(gain)) and np.all(np.isfinite(power_weight)) and np.all(np.isfinite(a_db))):
        raise ValueError("A-weighting evaluation exceeded float64 range")

    if is_scalar:
        return {
            "r_a": float(ra_f[0]),
            "r_a_1000": ra_1000,
            "linear_gain": float(gain[0]),
            "power_weight": float(power_weight[0]),
            "a_db": float(a_db[0]),
        }
    return {
        "r_a": ra_f,
        "r_a_1000": ra_1000,
        "linear_gain": gain,
        "power_weight": power_weight,
        "a_db": a_db,
    }


def differential_1f_noise_budget(
    frequencies_hz: np.ndarray,
    spacing_m: float,
    *,
    corner_frequency_hz: float = 1000.0,
    high_freq_floor_ratio: float = 0.5,
    sound_speed: float = 343.0,
) -> dict[str, object]:
    """Compare flat-floor Mic A vs 1/f Mic B under A-weighting, unweighted PSD,
    and two-mic endfire first-order differential equalization 1 / (2 * sin^2(2*pi*f*d/c)).

    Rejects spatial-aliasing zeros where sin(2*pi*f*d/c) == 0 on the evaluation grid.
    """
    freq = _finite_real_array(frequencies_hz, "frequencies_hz")
    if freq.ndim != 1 or freq.size < 1 or np.any(freq <= 0.0):
        raise ValueError("frequencies_hz must be a non-empty 1-D array of strictly positive frequencies")
    d_m = _finite_real_scalar(spacing_m, "spacing_m")
    f_c = _finite_real_scalar(corner_frequency_hz, "corner_frequency_hz")
    floor_ratio = _finite_real_scalar(high_freq_floor_ratio, "high_freq_floor_ratio")
    c_m_s = _finite_real_scalar(sound_speed, "sound_speed")
    if d_m <= 0.0 or f_c <= 0.0 or floor_ratio <= 0.0 or c_m_s <= 0.0:
        raise ValueError("spacing_m, corner_frequency_hz, high_freq_floor_ratio, and sound_speed must be positive")

    phase = 2.0 * np.pi * freq * d_m / c_m_s
    sin_term = np.sin(phase)
    if np.any(np.abs(sin_term) <= 1e-12):
        raise ValueError("differential array hits a spatial-aliasing zero sin(2*pi*f*d/c) == 0 on grid")

    psd_a = np.ones_like(freq, dtype=np.float64)
    psd_b = floor_ratio * (1.0 + f_c / freq)
    aw = iec61672_a_weighting(freq)
    w_a = np.asarray(aw["power_weight"], dtype=np.float64)

    diff_wng = 2.0 * sin_term * sin_term
    diff_gain_power = 1.0 / diff_wng

    unweighted_sum_a = float(np.sum(psd_a))
    unweighted_sum_b = float(np.sum(psd_b))
    a_weighted_sum_a = float(np.sum(psd_a * w_a))
    a_weighted_sum_b = float(np.sum(psd_b * w_a))
    diff_sum_a = float(np.sum(psd_a * diff_gain_power))
    diff_sum_b = float(np.sum(psd_b * diff_gain_power))

    f_min = float(np.min(freq))
    f_max = float(np.max(freq))
    asymp_int_a = (1.0 / f_min) - (1.0 / f_max)
    asymp_int_b = floor_ratio * (
        asymp_int_a + 0.5 * f_c * (1.0 / (f_min * f_min) - 1.0 / (f_max * f_max))
    )
    asymp_ratio = asymp_int_b / asymp_int_a if asymp_int_a > 0.0 else None

    return {
        "unweighted_ratio_b_over_a": unweighted_sum_b / unweighted_sum_a,
        "unweighted_ratio_db": float(10.0 * np.log10(unweighted_sum_b / unweighted_sum_a)),
        "a_weighted_ratio_b_over_a": a_weighted_sum_b / a_weighted_sum_a,
        "a_weighted_apparent_snr_advantage_b_db": float(-10.0 * np.log10(a_weighted_sum_b / a_weighted_sum_a)),
        "differential_grid_ratio_b_over_a": diff_sum_b / diff_sum_a,
        "differential_grid_penalty_b_db": float(10.0 * np.log10(diff_sum_b / diff_sum_a)),
        "differential_asymptotic_integral_ratio_b_over_a": (
            float(asymp_ratio) if asymp_ratio is not None else None
        ),
        "differential_asymptotic_penalty_b_db": (
            float(10.0 * np.log10(asymp_ratio)) if (asymp_ratio is not None and asymp_ratio > 0.0) else None
        ),
        "psd_b_over_a": psd_b,
        "differential_wng": diff_wng,
    }


def sofa_spherical_to_cartesian_and_colatitude(
    azimuth_sofa_deg: float,
    elevation_sofa_deg: float,
    radius_m: float = 1.0,
) -> dict[str, object]:
    """Convert AES69 SOFA spherical (az_ccw_deg, el_deg in [-90, 90], r_m > 0)
    to SOFA Cartesian (front, left, up), book angles/Cartesian, and colatitude [0, 180].

    Explicitly raises ValueError if ``elevation_sofa_deg`` is outside [-90.0, 90.0],
    ``radius_m <= 0``, or inputs are bool/complex/non-finite.
    """
    az_sofa = _finite_real_scalar(azimuth_sofa_deg, "azimuth_sofa_deg")
    el_sofa = _finite_real_scalar(elevation_sofa_deg, "elevation_sofa_deg")
    r_m = _finite_real_scalar(radius_m, "radius_m")
    if not (-90.0 <= el_sofa <= 90.0):
        raise ValueError("elevation_sofa_deg must lie within [-90.0, 90.0]")
    if r_m <= 0.0:
        raise ValueError("radius_m must be strictly positive")

    az_rad = math.radians(az_sofa)
    el_rad = math.radians(el_sofa)
    cos_el = math.cos(el_rad)
    x_sofa = r_m * cos_el * math.cos(az_rad)
    y_sofa = r_m * cos_el * math.sin(az_rad)
    z_sofa = r_m * math.sin(el_rad)

    colatitude_deg = 90.0 - el_sofa
    # Book convention: azimuth clockwise from +y (front) toward +x (right), wrapped to (-180, 180].
    book_az_deg = ((-az_sofa + 180.0) % 360.0) - 180.0
    if book_az_deg == -180.0 and -az_sofa > 0.0:
        book_az_deg = 180.0
    x_book = -y_sofa
    y_book = x_sofa
    z_book = z_sofa

    return {
        "azimuth_sofa_deg": az_sofa,
        "azimuth_sofa_wrapped_360_deg": az_sofa % 360.0,
        "elevation_sofa_deg": el_sofa,
        "radius_m": r_m,
        "colatitude_deg": colatitude_deg,
        "book_azimuth_deg": book_az_deg,
        "book_elevation_deg": el_sofa,
        "cartesian_sofa_m": [x_sofa, y_sofa, z_sofa],
        "cartesian_book_m": [x_book, y_book, z_book],
    }


def spatial_coherence_coloring_and_boundary() -> dict[str, object]:
    """E01-11: EVD vs Cholesky two-mic coloring, complex conjugate check, and 4-mic singular boundary."""
    rho_two = 24.0 / 25.0
    gamma_2 = equicorrelated_coherence(2, rho_two)
    evd_2 = hermitian_psd_sqrt_coloring(gamma_2)
    sqrt_2 = np.asarray(evd_2["coloring_matrix"], dtype=np.float64)
    chol_2 = np.linalg.cholesky(gamma_2)

    z_orth = np.array([[1.0, 1.0, -1.0, -1.0], [1.0, -1.0, 1.0, -1.0]], dtype=np.float64)
    v_evd = color_white_samples(sqrt_2, z_orth)
    v_chol = color_white_samples(chol_2, z_orth)
    cov_evd = (v_evd @ v_evd.T) / 4.0
    cov_chol = (v_chol @ v_chol.T) / 4.0

    gamma_2_complex = np.array([[1.0 + 0.0j, 0.96j], [-0.96j, 1.0 + 0.0j]], dtype=np.complex128)
    evd_2c = hermitian_psd_sqrt_coloring(gamma_2_complex)
    sqrt_2c = np.asarray(evd_2c["coloring_matrix"], dtype=np.complex128)
    hermitian_product = sqrt_2c @ sqrt_2c.conj().T
    wrong_transpose_product = sqrt_2c @ sqrt_2c.T

    rho_singular = -1.0 / 3.0
    gamma_4_singular = equicorrelated_coherence(4, rho_singular)
    evd_4_sing = hermitian_psd_sqrt_coloring(gamma_4_singular)
    sqrt_4_sing = np.asarray(evd_4_sing["coloring_matrix"], dtype=np.float64)
    z_4 = np.eye(4, dtype=np.float64)
    v_4_sing = color_white_samples(sqrt_4_sing, z_4)
    channel_sums = np.sum(v_4_sing, axis=0)

    # Exact sequential Cholesky pivots for 4-mic rho = -1/3: 4th pivot radicand vanishes identically.
    l11 = 1.0
    l21 = rho_singular
    l22 = math.sqrt(1.0 - l21 * l21)
    l31 = rho_singular
    l32 = (rho_singular - l31 * l21) / l22
    l33 = math.sqrt(1.0 - l31 * l31 - l32 * l32)
    l41 = rho_singular
    l42 = (rho_singular - l41 * l21) / l22
    l43 = (rho_singular - l41 * l31 - l42 * l32) / l33
    fourth_radicand = 1.0 - (l41 * l41 + l42 * l42 + l43 * l43)

    gamma_4_invalid = equicorrelated_coherence(4, -0.4)
    try:
        hermitian_psd_sqrt_coloring(gamma_4_invalid)
    except np.linalg.LinAlgError as err:
        invalid_rejected, invalid_reason = True, str(err)
    else:
        raise AssertionError("Expected LinAlgError for non-PSD rho = -0.4")

    return {
        "two_mic_real": {
            "rho": rho_two,
            "eigenvalues": np.asarray(evd_2["eigenvalues"]).tolist(),
            "evd_sqrt_matrix": sqrt_2.tolist(),
            "cholesky_lower_matrix": chol_2.tolist(),
            "evd_colored_samples": v_evd.tolist(),
            "cholesky_colored_samples": v_chol.tolist(),
            "evd_sample_covariance": cov_evd.tolist(),
            "cholesky_sample_covariance": cov_chol.tolist(),
        },
        "two_mic_quadrature_complex": {
            "sqrt_real": sqrt_2c.real.tolist(),
            "sqrt_imag": sqrt_2c.imag.tolist(),
            "hermitian_diagonal_real": np.diag(hermitian_product).real.tolist(),
            "unconjugated_transpose_diagonal_real": np.diag(wrong_transpose_product).real.tolist(),
        },
        "four_mic_singular_boundary": {
            "rho": rho_singular,
            "eigenvalues": np.asarray(evd_4_sing["eigenvalues"]).tolist(),
            "numerical_rank": int(evd_4_sing["numerical_rank"]),
            "evd_sqrt_matrix": sqrt_4_sing.tolist(),
            "cholesky_pivots_1_to_3": [l11, l22, l33],
            "cholesky_fourth_radicand": fourth_radicand,
            "colored_basis_column_sums": channel_sums.tolist(),
        },
        "non_psd_rejection_control": {
            "rho": -0.4,
            "min_eigenvalue_analytic": 1.0 + 3.0 * (-0.4),
            "rejected": invalid_rejected,
            "reason": invalid_reason,
        },
        "limits": (
            "Known target coherence matrix coloring at a single frequency bin; not blind noise covariance "
            "estimation or time-domain broadband FIR synthesis."
        ),
    }


def sensor_dba_trap_and_sofa_coordinates() -> dict[str, object]:
    """E01-12: IEC 61672-1 A-weighting vs 1/f differential noise budget and AES69 SOFA coordinates."""
    spot_freqs = np.array([50.0, 100.0, 200.0, 500.0, 1000.0, 2500.0], dtype=np.float64)
    spot_aw = iec61672_a_weighting(spot_freqs)
    spot_rows = [
        {
            "frequency_hz": float(f),
            "a_db": float(db),
            "power_weight": float(pw),
        }
        for f, db, pw in zip(
            spot_freqs,
            np.asarray(spot_aw["a_db"]),
            np.asarray(spot_aw["power_weight"]),
        )
    ]

    grid_dba = np.arange(50.0, 8000.0 + 1.0, 10.0, dtype=np.float64)
    budget_dba_band = differential_1f_noise_budget(grid_dba, 0.01)

    grid_speech = np.arange(100.0, 4000.0 + 1.0, 10.0, dtype=np.float64)
    budget_speech_band = differential_1f_noise_budget(grid_speech, 0.01)

    spot_100_1000 = differential_1f_noise_budget(np.array([100.0, 1000.0], dtype=np.float64), 0.01)
    wng_100 = float(np.asarray(spot_100_1000["differential_wng"])[0])
    wng_1000 = float(np.asarray(spot_100_1000["differential_wng"])[1])

    sofa_ref = sofa_spherical_to_cartesian_and_colatitude(-30.0, 30.0, 1.0)
    # Misuse diagnostics: passing elevation delta directly into spherical-harmonic colatitude vartheta.
    misuse_horizontal_sin_theta = math.sin(math.radians(0.0))
    misuse_lower_hemisphere_sin_theta = math.sin(math.radians(-30.0))
    misuse_lower_hemisphere_cos_theta = math.cos(math.radians(-30.0))

    return {
        "iec61672_spot_table": spot_rows,
        "r_a_1000": float(spot_aw["r_a_1000"]),
        "mic_b_over_a_unweighted_100hz": float(0.5 * (1.0 + 1000.0 / 100.0)),
        "mic_b_over_a_unweighted_100hz_db": float(10.0 * math.log10(5.5)),
        "dba_band_50_8000_hz": {
            "num_points": int(grid_dba.size),
            "a_weighted_ratio_b_over_a": float(budget_dba_band["a_weighted_ratio_b_over_a"]),
            "apparent_dba_advantage_b_db": float(budget_dba_band["a_weighted_apparent_snr_advantage_b_db"]),
        },
        "endfire_differential_1cm": {
            "wng_100hz": wng_100,
            "wng_100hz_db": float(10.0 * math.log10(wng_100)),
            "wng_1000hz": wng_1000,
            "wng_1000hz_db": float(10.0 * math.log10(wng_1000)),
            "speech_band_100_4000_grid_ratio_b_over_a": float(
                budget_speech_band["differential_grid_ratio_b_over_a"]
            ),
            "speech_band_100_4000_grid_penalty_b_db": float(
                budget_speech_band["differential_grid_penalty_b_db"]
            ),
            "speech_band_100_4000_asymptotic_ratio_b_over_a": float(
                budget_speech_band["differential_asymptotic_integral_ratio_b_over_a"]
            ),
            "speech_band_100_4000_asymptotic_penalty_b_db": float(
                budget_speech_band["differential_asymptotic_penalty_b_db"]
            ),
        },
        "sofa_reference_point": sofa_ref,
        "colatitude_misuse_diagnostics": {
            "horizontal_delta_0_sin": misuse_horizontal_sin_theta,
            "lower_delta_minus30_sin": misuse_lower_hemisphere_sin_theta,
            "lower_delta_minus30_z_sign": misuse_lower_hemisphere_cos_theta,
        },
        "limits": (
            "Analytical parametric self-noise PSD model and AES69 SOFA coordinate definitions; "
            "not anechoic chamber MEMS characterization or measured HRIR dataset parsing."
        ),
    }
