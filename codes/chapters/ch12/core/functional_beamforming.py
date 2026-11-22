"""Dougherty functional beamforming (R^(1/nu)) teaching kernels for Ch 12.

Positive-frequency convention matches ``codes.chapters.ch12.core.imaging``:
``A`` has shape ``(M, G)`` (microphones x grid), ``R`` is an ``(M, M)``
Hermitian positive semidefinite cross-spectral matrix (full CSM, never
diagonal-removed), and Formulation III ('true level' in locked Acoular) uses:

    u_i = a_i / ||a_i||_2,
    b_FB,i(nu) = (u_i^H R^(1/nu) u_i)^nu / ||a_i||_2^2,   nu >= 1.

At nu = 1 this reduces identically to conventional true-level beamforming
w_i^H R w_i with w_i = a_i / ||a_i||_2^2. No SciPy or Acoular import is used.
"""
from __future__ import annotations

import math
import numpy as np


def _validate_real_scalar(value: object, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(
        value, (int, float, np.integer, np.floating)
    ):
        raise ValueError(f"{name} must be a real numeric scalar, not bool")
    scalar = float(value)
    if not math.isfinite(scalar):
        raise ValueError(f"{name} must be finite")
    return scalar


def _validate_complex_matrix(value: object, name: str) -> np.ndarray:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be numeric, not bool")
    raw = np.asarray(value)
    if raw.dtype.kind not in "iufc":
        raise ValueError(f"{name} must be numeric, not bool or text")
    if isinstance(value, (list, tuple)):
        def _check_no_bool(items: object) -> None:
            if isinstance(items, (list, tuple)):
                for item in items:
                    _check_no_bool(item)
            elif isinstance(items, (bool, np.bool_)):
                raise ValueError(f"{name} must not contain bool")

        _check_no_bool(value)
    with np.errstate(over="ignore", invalid="ignore"):
        array = np.asarray(raw, dtype=np.complex128)
    if array.ndim != 2 or min(array.shape) < 1:
        raise ValueError(f"{name} must be a non-empty 2D matrix")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


def _validate_hpsd_csm(value: object, *, tol: float = 1e-12) -> np.ndarray:
    tol_val = _validate_real_scalar(tol, "tol")
    if tol_val <= 0.0:
        raise ValueError("tol must be positive")
    array = _validate_complex_matrix(value, "csm")
    if array.shape[0] != array.shape[1]:
        raise ValueError("csm must be square")
    scale = float(np.max(np.abs(array)))
    if scale == 0.0:
        return np.zeros_like(array, dtype=np.complex128)
    normalized = array.real / scale + 1j * (array.imag / scale)
    if float(np.max(np.abs(normalized - normalized.conj().T))) > tol_val:
        raise ValueError("csm must be Hermitian")
    sym_norm = 0.5 * (normalized + normalized.conj().T)
    min_eig = float(np.min(np.linalg.eigvalsh(sym_norm)))
    if min_eig < -tol_val:
        raise ValueError("csm must be positive semidefinite (full CSM required)")
    if np.array_equal(array, array.conj().T):
        return array.copy()
    return 0.5 * (array + array.conj().T)


def hermitian_fractional_power(
    csm: object,
    exponent: float,
    *,
    tol: float = 1e-12,
) -> np.ndarray:
    """Compute H^exponent for Hermitian positive semidefinite H and exponent > 0."""
    exp_val = _validate_real_scalar(exponent, "exponent")
    if exp_val <= 0.0:
        raise ValueError("exponent must be strictly positive")
    herm = _validate_hpsd_csm(csm, tol=tol)
    scale = float(np.max(np.abs(herm)))
    if scale == 0.0:
        return np.zeros_like(herm, dtype=np.complex128)
    herm_unit = herm / scale
    eigenvalues, eigenvectors = np.linalg.eigh(herm_unit)
    eig_max = max(float(np.max(np.abs(eigenvalues))), 1e-300)
    if float(np.min(eigenvalues)) < -tol * eig_max:
        raise ValueError("csm must be positive semidefinite")
    clipped = np.where(eigenvalues > tol * eig_max, eigenvalues, 0.0)
    with np.errstate(over="raise", invalid="raise"):
        powered = (scale ** exp_val) * np.power(clipped, exp_val)
        reconstructed = (eigenvectors * powered) @ eigenvectors.conj().T
    if not np.all(np.isfinite(reconstructed)):
        raise ValueError("fractional matrix power exceeds float64 range")
    return 0.5 * (reconstructed + reconstructed.conj().T)


def closed_form_2x2_sqrt(csm: object, *, tol: float = 1e-12) -> np.ndarray:
    """Exact 2x2 Cayley-Hamilton Hermitian positive semidefinite square root.

    For 2x2 H >= 0:
        H^(1/2) = (H + sqrt(det H) I_2) / sqrt(tr H + 2 sqrt(det H)).
    """
    herm = _validate_hpsd_csm(csm, tol=tol)
    if herm.shape != (2, 2):
        raise ValueError("closed_form_2x2_sqrt requires a 2x2 matrix")
    scale = float(np.max(np.abs(herm)))
    if scale == 0.0:
        return np.zeros((2, 2), dtype=np.complex128)
    herm_unit = herm / scale
    trace_val = float(herm_unit[0, 0].real + herm_unit[1, 1].real)
    if trace_val <= 0.0:
        return np.zeros((2, 2), dtype=np.complex128)
    off_sq = float(herm_unit[0, 1].real ** 2 + herm_unit[0, 1].imag ** 2)
    det_val = float(herm_unit[0, 0].real * herm_unit[1, 1].real - off_sq)
    if det_val < -tol * max(trace_val * trace_val, 1e-300):
        raise ValueError("2x2 csm determinant is negative beyond tolerance")
    sqrt_det = math.sqrt(max(det_val, 0.0))
    denom = math.sqrt(max(trace_val + 2.0 * sqrt_det, 0.0))
    if denom <= 0.0:
        return np.zeros((2, 2), dtype=np.complex128)
    root = math.sqrt(scale) * (herm_unit + sqrt_det * np.eye(2, dtype=np.complex128)) / denom
    return 0.5 * (root + root.conj().T)


def functional_beamforming_map(
    csm: object,
    steering: object,
    nu: float = 1.0,
    *,
    tol: float = 1e-12,
) -> np.ndarray:
    """Evaluate Formulation-III ('true level') functional beamforming map.

    Parameters
    ----------
    csm : (M, M) array_like
        Full Hermitian positive semidefinite cross-spectral matrix.
    steering : (M, G) array_like
        Propagation/steering vectors a_i for G scan grid points.
    nu : float, default 1.0
        Order parameter nu >= 1.0 (``gamma`` in Acoular ``BeamformerFunctional``).
    """
    nu_val = _validate_real_scalar(nu, "nu")
    if nu_val < 1.0:
        raise ValueError("functional beamforming requires nu >= 1.0")
    herm = _validate_hpsd_csm(csm, tol=tol)
    a_mat = _validate_complex_matrix(steering, "steering")
    if a_mat.shape[0] != herm.shape[0]:
        raise ValueError("steering row count must match csm dimension")
    norms_sq = np.sum(np.abs(a_mat) ** 2, axis=0)
    if np.any(norms_sq <= 0.0) or not np.all(np.isfinite(norms_sq)):
        raise ValueError("steering columns must have positive finite squared norm")
    csm_root = hermitian_fractional_power(herm, 1.0 / nu_val, tol=tol)
    unit_steering = a_mat / np.sqrt(norms_sq)
    quad = np.einsum("mg,mn,ng->g", unit_steering.conj(), csm_root, unit_steering).real
    quad_scale = max(1.0, float(np.max(np.abs(quad))))
    if float(np.min(quad)) < -tol * quad_scale:
        raise ValueError("fractional CSM quadratic form became negative")
    clipped_quad = np.maximum(quad, 0.0)
    with np.errstate(over="raise", invalid="raise"):
        out = np.power(clipped_quad, nu_val) / norms_sq
    if not np.all(np.isfinite(out)):
        raise ValueError("functional beamforming output exceeds float64 range")
    return np.asarray(out, dtype=np.float64)


def functional_beamforming_controls() -> dict[str, dict[str, object]]:
    """Return deterministic E12-17 and E12-18 teaching controls."""
    # Shared two-cell benchmark from Chapter 12:
    # a_1 = [1, 1]^T, a_2 = [1, exp(-j 2 pi / 3)]^T, P_11 = P_22 = 1, P_12 = P_21 = 1/4.
    a_two = np.array(
        [[1.0, 1.0], [1.0, np.exp(-1j * 2.0 * np.pi / 3.0)]],
        dtype=np.complex128,
    )
    psf_21 = 0.25

    # --- E12-17: Single-source peak invariance & exponential sidelobe suppression ---
    q_single = 1.0
    r_single = q_single * np.outer(a_two[:, 0], a_two[:, 0].conj())
    single_cases = []
    for nu in (1, 2, 4):
        fb_vals = functional_beamforming_map(r_single, a_two, nu=float(nu))
        analytic_vals = np.array([q_single, q_single * (psf_21 ** nu)], dtype=np.float64)
        sidelobe_db = float(10.0 * math.log10(fb_vals[1] / fb_vals[0]))
        single_cases.append(
            {
                "nu": nu,
                "fb_map": fb_vals,
                "analytic_map": analytic_vals,
                "max_abs_error": float(np.max(np.abs(fb_vals - analytic_vals))),
                "sidelobe_db": sidelobe_db,
            }
        )

    # Near-field unequal-norm single-source check: ||a_1||_2^2 = 5, ||a_2||_2^2 = 1.5
    a_unequal = np.array([[1.0, 1.0], [2.0, -0.5 + 0.5j]], dtype=np.complex128)
    q_near = 1.5
    r_near = q_near * np.outer(a_unequal[:, 0], a_unequal[:, 0].conj())
    norm1_sq = float(np.vdot(a_unequal[:, 0], a_unequal[:, 0]).real)
    norm2_sq = float(np.vdot(a_unequal[:, 1], a_unequal[:, 1]).real)
    cos_sq = float(
        abs(np.vdot(a_unequal[:, 1], a_unequal[:, 0])) ** 2 / (norm1_sq * norm2_sq)
    )
    unequal_cases = []
    for nu in (1, 2, 4):
        fb_near = functional_beamforming_map(r_near, a_unequal, nu=float(nu))
        analytic_near = np.array(
            [q_near, q_near * (norm1_sq / norm2_sq) * (cos_sq ** nu)],
            dtype=np.float64,
        )
        unequal_cases.append(
            {
                "nu": nu,
                "fb_map": fb_near,
                "analytic_map": analytic_near,
                "peak_error": float(abs(fb_near[0] - q_near)),
                "max_abs_error": float(np.max(np.abs(fb_near - analytic_near))),
            }
        )

    removed_diag_r = r_single - np.diag(np.diag(r_single))
    removed_eigs = np.linalg.eigvalsh(removed_diag_r)
    diag_removal_rejected = False
    try:
        functional_beamforming_map(removed_diag_r, a_two, nu=2.0)
    except ValueError:
        diag_removal_rejected = True
    if not diag_removal_rejected:
        raise AssertionError("expected diagonal-removed indefinite CSM to be rejected")

    e12_17 = {
        "A_real_imag": a_two,
        "single_source_q": [1.0, 0.0],
        "single_source_R_real_imag": r_single,
        "nu_cases": single_cases,
        "unequal_norm_A_real_imag": a_unequal,
        "unequal_norm_q1": q_near,
        "unequal_norm_cases": unequal_cases,
        "diagonal_removed_eigenvalues": removed_eigs,
        "diagonal_removed_rejected": diag_removal_rejected,
        "scope": (
            "Formulation-III true-level functional beamforming on full HPSD CSM; "
            "single-source peak is invariant for all nu >= 1 while sidelobes scale as P_i1^nu"
        ),
    }

    # --- E12-18: 2x2 closed-form sqrt, power non-conservation & off-grid weak source ---
    q_two = np.array([1.0, 0.25], dtype=np.float64)
    r_two = a_two @ np.diag(q_two) @ a_two.conj().T
    sqrt_closed = closed_form_2x2_sqrt(r_two)
    sqrt_eigh = hermitian_fractional_power(r_two, 0.5)
    b_conv = functional_beamforming_map(r_two, a_two, nu=1.0)
    b_fb2 = functional_beamforming_map(r_two, a_two, nu=2.0)

    sqrt3 = math.sqrt(3.0)
    b_fb2_analytic = np.array(
        [(869.0 + 6.0 * sqrt3) / 832.0, (11.0 + 6.0 * sqrt3) / 52.0],
        dtype=np.float64,
    )
    psf_sum_nu1 = 1.25  # 1 + 1/4
    psf_sum_nu2 = 17.0 / 16.0  # 1 + 1/16
    fb2_raw_sum = float(np.sum(b_fb2))
    fb2_raw_sum_analytic = (1045.0 + 102.0 * sqrt3) / 832.0
    fb2_psf_norm_sum = fb2_raw_sum / psf_sum_nu2
    fb2_psf_norm_sum_analytic = (1045.0 + 102.0 * sqrt3) / 884.0
    hadamard_sq_psf_linear_map = np.array([[1.0, 1.0 / 16.0], [1.0 / 16.0, 1.0]]) @ q_two

    # M = 3 off-grid / steering-mismatch weak-source underestimation control:
    # a_1 = [1, 1, 1]^T (q_1 = 1), a_true,2 = [1, z, z^2]^T (q_2 = 0.25),
    # u_3 = [1, z^*, (z^*)^2]^T / sqrt(3) in ker(R), eta = |u_scan,2^H u_true,2|^2 = 4/5.
    z_phase = np.exp(-1j * 2.0 * np.pi / 3.0)
    a3_strong = np.ones(3, dtype=np.complex128)
    a3_weak_true = np.array([1.0, z_phase, z_phase ** 2], dtype=np.complex128)
    u3_null = np.array([1.0, z_phase.conj(), (z_phase.conj()) ** 2], dtype=np.complex128) / math.sqrt(3.0)
    eta = 0.8
    a3_weak_scan = math.sqrt(3.0) * (
        math.sqrt(eta) * (a3_weak_true / math.sqrt(3.0))
        + math.sqrt(1.0 - eta) * u3_null
    )
    a3_scan_grid = np.column_stack([a3_strong, a3_weak_scan])
    r3_two = (
        1.0 * np.outer(a3_strong, a3_strong.conj())
        + 0.25 * np.outer(a3_weak_true, a3_weak_true.conj())
    )
    offgrid_cases = []
    for nu in (1, 2, 8):
        vals = functional_beamforming_map(r3_two, a3_scan_grid, nu=float(nu))
        analytic_weak = 0.25 * (eta ** nu)
        offgrid_cases.append(
            {
                "nu": nu,
                "scan_map": vals,
                "expected_strong": 1.0,
                "expected_weak": analytic_weak,
                "weak_error_db_re_q2": float(10.0 * math.log10(vals[1] / 0.25)),
                "max_abs_error": float(
                    max(abs(vals[0] - 1.0), abs(vals[1] - analytic_weak))
                ),
            }
        )

    # Destructive coherence control: x = a_1 - 0.5 * exp(j*pi/3) * a_2 (|s_2|^2 = 0.25, a_2^H x = 0) => b_FB,2(nu) = 0 < q_2 = 0.25
    x_coherent = a_two @ np.array([1.0, -0.5 * np.exp(1j * np.pi / 3.0)], dtype=np.complex128)
    r_coherent = np.outer(x_coherent, x_coherent.conj())
    coherent_fb2 = functional_beamforming_map(r_coherent, a_two, nu=2.0)

    e12_18 = {
        "q_true": q_two,
        "R_real_imag": r_two,
        "trace": float(np.trace(r_two).real),
        "determinant": float(np.linalg.det(r_two).real),
        "eigenvalues": np.linalg.eigvalsh(r_two),
        "sqrt_closed_form_real_imag": sqrt_closed,
        "sqrt_eigh_real_imag": sqrt_eigh,
        "sqrt_method_max_abs_error": float(np.max(np.abs(sqrt_closed - sqrt_eigh))),
        "sqrt_reconstruction_max_abs_error": float(
            np.max(np.abs(sqrt_closed @ sqrt_closed - r_two))
        ),
        "conventional_map_nu_1": b_conv,
        "fb_map_nu_2": b_fb2,
        "fb_map_nu_2_analytic": b_fb2_analytic,
        "fb_map_nu_2_max_abs_error": float(np.max(np.abs(b_fb2 - b_fb2_analytic))),
        "power_accounting": {
            "true_source_sum": float(np.sum(q_two)),
            "conventional_raw_sum": float(np.sum(b_conv)),
            "conventional_psf_normalized_sum": float(np.sum(b_conv) / psf_sum_nu1),
            "fb_nu_2_raw_sum": fb2_raw_sum,
            "fb_nu_2_raw_sum_analytic": fb2_raw_sum_analytic,
            "fb_nu_2_psf_column_sum": psf_sum_nu2,
            "fb_nu_2_psf_normalized_sum": fb2_psf_norm_sum,
            "fb_nu_2_psf_normalized_sum_analytic": fb2_psf_norm_sum_analytic,
            "hadamard_squared_psf_linear_map": hadamard_sq_psf_linear_map,
            "hadamard_squared_psf_linear_sum": float(np.sum(hadamard_sq_psf_linear_map)),
        },
        "offgrid_weak_source_3mic": {
            "eta_squared_cosine": eta,
            "csm_eigenvalues": np.linalg.eigvalsh(r3_two),
            "nu_cases": offgrid_cases,
        },
        "coherent_dipole_fb_nu_2": coherent_fb2,
        "scope": (
            "matrix fractional powers are non-additive: FB does not conserve integrated "
            "multi-source power and suppresses off-grid weak sources as eta^nu when rank(R) < M"
        ),
    }

    return {"E12-17": e12_17, "E12-18": e12_18}
