"""Pure-NumPy teaching kernels for Hammerstein bases and stereo/loop AEC controls.

Covers Chapter 06 exercises E06-29 (odd power vs orthogonal Legendre Hammerstein
conditioning and coordinate decoupling on U(-1, 1)) and E06-30 (LTI vs
time-varying first-order all-pass stereo decorrelation, physical ANC/ITU-T
P.1150 ICC delay deadlines vs digital AEC buffering, and foreground/shadow
dual-filter switching). No SciPy or external dependencies; all computations are
deterministic float64.
"""
from __future__ import annotations

import math
import numbers
import numpy as np


def _finite_real_array(value, name: str = "input") -> np.ndarray:
    original = np.asarray(value)
    if original.dtype.kind not in "iuf":
        raise ValueError(f"{name} must contain real numeric values, not complex, bool or text")
    with np.errstate(over="ignore", invalid="ignore"):
        array = np.asarray(original, dtype=np.float64)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite and representable as float64")
    return array


def _finite_real_scalar(value, name: str = "input") -> float:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a real scalar, not bool")
    array = _finite_real_array(value, name)
    if array.ndim != 0:
        raise ValueError(f"{name} must be a scalar")
    return float(array)


def _integer(value, name: str, *, minimum: int | None = None) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, numbers.Integral):
        raise ValueError(f"{name} must be an integer")
    result = int(value)
    if minimum is not None and result < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return result


def odd_power_basis(x: np.ndarray) -> np.ndarray:
    """Evaluate 3-branch odd power basis ``[x, x^3, x^5]`` of shape ``(N, 3)``."""
    samples = _finite_real_array(x, "x")
    if samples.ndim != 1 or samples.size == 0:
        raise ValueError("x must be a non-empty 1-D real array")
    with np.errstate(over="raise", invalid="raise"):
        try:
            x2 = samples * samples
            x3 = x2 * samples
            x5 = x3 * x2
        except FloatingPointError as error:
            raise ValueError("odd_power_basis overflows float64") from error
    basis = np.column_stack((samples, x3, x5))
    if not np.all(np.isfinite(basis)):
        raise ValueError("odd_power_basis produced non-finite values")
    return basis


def odd_legendre_basis(x: np.ndarray) -> np.ndarray:
    """Evaluate 3-branch odd Legendre basis ``[P_1(x), P_3(x), P_5(x)]`` of shape ``(N, 3)``."""
    samples = _finite_real_array(x, "x")
    if samples.ndim != 1 or samples.size == 0:
        raise ValueError("x must be a non-empty 1-D real array")
    with np.errstate(over="raise", invalid="raise"):
        try:
            x2 = samples * samples
            x3 = x2 * samples
            x5 = x3 * x2
            p1 = samples
            p3 = 0.5 * (5.0 * x3 - 3.0 * samples)
            p5 = 0.125 * (63.0 * x5 - 70.0 * x3 + 15.0 * samples)
        except FloatingPointError as error:
            raise ValueError("odd_legendre_basis overflows float64") from error
    basis = np.column_stack((p1, p3, p5))
    if not np.all(np.isfinite(basis)):
        raise ValueError("odd_legendre_basis produced non-finite values")
    return basis


def power_to_legendre_weights(power_coefficients: np.ndarray) -> np.ndarray:
    """Map ``[c_1, c_3, c_5]`` in ``c_1 x + c_3 x^3 + c_5 x^5`` to ``[w_1, w_3, w_5]`` on ``[P_1, P_3, P_5]``."""
    coeffs = _finite_real_array(power_coefficients, "power_coefficients")
    if coeffs.shape != (3,):
        raise ValueError("power_coefficients must have shape (3,)")
    transform = np.array(
        [
            [1.0, 3.0 / 5.0, 3.0 / 7.0],
            [0.0, 2.0 / 5.0, 4.0 / 9.0],
            [0.0, 0.0, 8.0 / 63.0],
        ],
        dtype=np.float64,
    )
    weights = transform @ coeffs
    if not np.all(np.isfinite(weights)):
        raise ValueError("transformed Legendre weights are non-finite")
    return weights


def hammerstein_orthogonal_basis() -> dict:
    """E06-29: Conditioning and coordinate decoupling of odd power vs Legendre bases on U(-1, 1)."""
    orders = np.array([1, 3, 5], dtype=np.int64)
    raw_power_gram = 1.0 / (orders[:, None] + orders[None, :] + 1.0)
    raw_eigvals = np.linalg.eigvalsh(raw_power_gram)
    raw_cond_3 = float(raw_eigvals[-1] / raw_eigvals[0])

    raw_power_2 = raw_power_gram[:2, :2]
    raw_eigvals_2 = np.linalg.eigvalsh(raw_power_2)
    raw_cond_2 = float(raw_eigvals_2[-1] / raw_eigvals_2[0])

    scale_diag = np.sqrt(np.array([3.0, 7.0, 11.0], dtype=np.float64))
    unit_power_gram = (scale_diag[:, None] * raw_power_gram) * scale_diag[None, :]
    unit_eigvals = np.linalg.eigvalsh(unit_power_gram)
    unit_cond_3 = float(unit_eigvals[-1] / unit_eigvals[0])

    unit_power_2 = unit_power_gram[:2, :2]
    unit_eigvals_2 = np.linalg.eigvalsh(unit_power_2)
    unit_cond_2 = float(unit_eigvals_2[-1] / unit_eigvals_2[0])

    legendre_gram = np.diag([1.0 / 3.0, 1.0 / 7.0, 1.0 / 11.0])
    legendre_eigvals = np.linalg.eigvalsh(legendre_gram)
    legendre_cond_3 = float(legendre_eigvals[-1] / legendre_eigvals[0])
    orthonormal_legendre_gram = np.eye(3, dtype=np.float64)

    true_power_coeffs = np.array([1.0, -0.5, 0.63], dtype=np.float64)
    basis_transform = np.array(
        [
            [1.0, 3.0 / 5.0, 3.0 / 7.0],
            [0.0, 2.0 / 5.0, 4.0 / 9.0],
            [0.0, 0.0, 8.0 / 63.0],
        ],
        dtype=np.float64,
    )
    legendre_weights = power_to_legendre_weights(true_power_coeffs)

    # Optimal L2([-1, 1]) single-branch linear fit E[x y(x)] / E[x^2].
    cross_power = raw_power_gram @ true_power_coeffs
    truncated_linear_gain = float(cross_power[0] / raw_power_gram[0, 0])
    total_output_variance = float(true_power_coeffs @ cross_power)
    linear_explained_variance = float((legendre_weights[0] ** 2) / 3.0)
    residual_nonlinear_variance = float(
        (legendre_weights[1] ** 2) / 7.0 + (legendre_weights[2] ** 2) / 11.0
    )

    grid = np.linspace(-1.0, 1.0, 9)
    y_power = odd_power_basis(grid) @ true_power_coeffs
    y_legendre = odd_legendre_basis(grid) @ legendre_weights
    max_pointwise_discrepancy = float(np.max(np.abs(y_power - y_legendre)))

    return {
        "excitation_distribution": "uniform_minus1_to_plus1",
        "odd_orders": orders.tolist(),
        "raw_power_gram": raw_power_gram.tolist(),
        "raw_power_determinant": float(np.linalg.det(raw_power_gram)),
        "raw_power_eigenvalues_ascending": raw_eigvals.tolist(),
        "raw_power_condition_number_3branch": raw_cond_3,
        "raw_power_condition_number_2branch": raw_cond_2,
        "unit_variance_scales": scale_diag.tolist(),
        "unit_power_gram": unit_power_gram.tolist(),
        "unit_power_off_diagonals": {
            "rho_13": float(unit_power_gram[0, 1]),
            "rho_15": float(unit_power_gram[0, 2]),
            "rho_35": float(unit_power_gram[1, 2]),
        },
        "unit_power_determinant": float(np.linalg.det(unit_power_gram)),
        "unit_power_eigenvalues_ascending": unit_eigvals.tolist(),
        "unit_power_condition_number_3branch": unit_cond_3,
        "unit_power_condition_number_2branch": unit_cond_2,
        "legendre_gram": legendre_gram.tolist(),
        "legendre_condition_number_3branch": legendre_cond_3,
        "orthonormal_legendre_gram": orthonormal_legendre_gram.tolist(),
        "orthonormal_legendre_condition_number_3branch": 1.0,
        "power_to_legendre_transform": basis_transform.tolist(),
        "true_power_coefficients": true_power_coeffs.tolist(),
        "legendre_coefficients": legendre_weights.tolist(),
        "best_linear_truncation_gain": truncated_linear_gain,
        "linear_gain_matches_first_legendre_weight": bool(
            abs(truncated_linear_gain - float(legendre_weights[0])) <= 1e-15
        ),
        "total_output_variance": total_output_variance,
        "linear_explained_variance": linear_explained_variance,
        "residual_nonlinear_variance": residual_nonlinear_variance,
        "grid_points": grid.tolist(),
        "max_pointwise_basis_difference": max_pointwise_discrepancy,
        "scope": (
            "memoryless odd polynomial distortion under uniform [-1, 1] amplitude; "
            "real speech has non-uniform amplitude and temporal correlation"
        ),
    }


def first_order_allpass_decorrelate(
    source: np.ndarray,
    pole_trajectory: np.ndarray,
) -> np.ndarray:
    """Causal 1st-order real all-pass ``x_2[n] = -a_n x_2[n-1] + a_n s[n] + s[n-1]``."""
    s = _finite_real_array(source, "source")
    a = _finite_real_array(pole_trajectory, "pole_trajectory")
    if s.ndim != 1 or s.size == 0 or a.shape != s.shape:
        raise ValueError("source and pole_trajectory must be equal-length non-empty 1-D arrays")
    if np.any(np.abs(a) >= 1.0):
        raise ValueError("pole_trajectory magnitudes must be strictly less than 1")
    out = np.zeros_like(s)
    prev_s = 0.0
    prev_y = 0.0
    for idx in range(s.size):
        an = float(a[idx])
        sn = float(s[idx])
        yn = -an * prev_y + an * sn + prev_s
        if not math.isfinite(yn):
            raise ValueError("all-pass recursion overflowed float64")
        out[idx] = yn
        prev_s = sn
        prev_y = yn
    return out


def two_tap_stereo_regressor_matrix(
    x1: np.ndarray,
    x2: np.ndarray,
    sample_indices: np.ndarray | list[int] | tuple[int, ...],
) -> np.ndarray:
    """Build snapshot rows ``[x1[n], x1[n-1], x2[n], x2[n-1]]`` with zero pre-history."""
    left = _finite_real_array(x1, "x1")
    right = _finite_real_array(x2, "x2")
    if left.ndim != 1 or left.size == 0 or right.shape != left.shape:
        raise ValueError("x1 and x2 must be equal-length non-empty 1-D arrays")
    if isinstance(sample_indices, (str, bytes)) or not hasattr(sample_indices, "__iter__"):
        raise ValueError("sample_indices must be an iterable of sample indices")
    indices = [_integer(idx, "sample_index", minimum=0) for idx in sample_indices]
    if not indices or max(indices) >= left.size:
        raise ValueError("sample_indices must be non-empty and within signal bounds")
    rows = np.empty((len(indices), 4), dtype=np.float64)
    for row_idx, n in enumerate(indices):
        rows[row_idx, 0] = left[n]
        rows[row_idx, 1] = left[n - 1] if n >= 1 else 0.0
        rows[row_idx, 2] = right[n]
        rows[row_idx, 3] = right[n - 1] if n >= 1 else 0.0
    return rows


def dual_filter_switch_decision(
    foreground_error_energy: float,
    shadow_error_energy: float,
    coherence: float,
    *,
    copy_ratio: float = 0.25,
    reset_ratio: float = 2.0,
    coherence_threshold: float = 0.8,
) -> dict:
    """Evaluate foreground/shadow AEC switching rule (Eq. 6-43)."""
    e_fg = _finite_real_scalar(foreground_error_energy, "foreground_error_energy")
    e_sh = _finite_real_scalar(shadow_error_energy, "shadow_error_energy")
    coh = _finite_real_scalar(coherence, "coherence")
    eta_copy = _finite_real_scalar(copy_ratio, "copy_ratio")
    eta_reset = _finite_real_scalar(reset_ratio, "reset_ratio")
    tau_c = _finite_real_scalar(coherence_threshold, "coherence_threshold")
    if e_fg <= 0.0 or e_sh < 0.0:
        raise ValueError("foreground_error_energy must be > 0 and shadow_error_energy >= 0")
    if not 0.0 <= coh <= 1.0:
        raise ValueError("coherence must lie in [0, 1]")
    if not (0.0 < eta_copy < 1.0 < eta_reset) or not 0.0 < tau_c <= 1.0:
        raise ValueError("require 0 < copy_ratio < 1 < reset_ratio and 0 < coherence_threshold <= 1")
    ratio = e_sh / e_fg
    if (e_sh < eta_copy * e_fg) and (coh >= tau_c):
        action = "copy_shadow_to_foreground"
    elif e_sh > eta_reset * e_fg:
        action = "reset_shadow_from_foreground"
    else:
        action = "hold"
    return {
        "foreground_error_energy": e_fg,
        "shadow_error_energy": e_sh,
        "shadow_to_foreground_ratio": ratio,
        "coherence": coh,
        "action": action,
    }


def stereo_decorrelation_and_loop_deadlines() -> dict:
    """E06-30: LTI vs TV-APF stereo identifiability, ANC/ICC deadlines, and dual-filter switching."""
    source = np.array([1.0, 0.0, 0.0, 1.0, 0.0], dtype=np.float64)
    snapshot_indices = [0, 1, 3, 4]
    true_stereo_path = np.array([0.8, 0.2, 0.4, 0.1], dtype=np.float64)

    lti_pole = 0.5
    x2_lti = first_order_allpass_decorrelate(source, np.full(source.size, lti_pole))
    u_lti = two_tap_stereo_regressor_matrix(source, x2_lti, snapshot_indices)
    lti_null_vector = np.array([-lti_pole, -1.0, 1.0, lti_pole], dtype=np.float64)
    lti_null_residual = u_lti @ lti_null_vector

    tv_poles = np.array([0.0, 0.0, 0.0, 0.5, 0.5], dtype=np.float64)
    x2_tv = first_order_allpass_decorrelate(source, tv_poles)
    u_tv = two_tap_stereo_regressor_matrix(source, x2_tv, snapshot_indices)
    y_tv = u_tv @ true_stereo_path
    recovered_stereo_path = np.linalg.solve(u_tv, y_tv)

    # Physical ANC vs digital open-loop AEC vs ITU-T P.1150 ICC timing budgets.
    c_mps = 343.0
    fs_hz = 16000.0
    d_primary_m = 0.12
    d_secondary_m = 0.02
    tau_acoust_ms = 1000.0 * d_primary_m / c_mps
    tau_sec_ms = 1000.0 * d_secondary_m / c_mps
    tau_elec_max_ms = tau_acoust_ms - tau_sec_ms
    tau_elec_max_samples = tau_elec_max_ms * 1e-3 * fs_hz
    pbfdaf_block_samples = 128
    pbfdaf_block_ms = 1000.0 * pbfdaf_block_samples / fs_hz
    aec_buffer_delay_samples = 800
    aec_buffer_delay_ms = 1000.0 * aec_buffer_delay_samples / fs_hz

    icc_alpha = 0.8
    icc_p1150_limit_ms = 15.0
    icc_cases = []
    for loop_delay_ms in (10.0, 25.0):
        loop_delay_s = loop_delay_ms * 1e-3
        notch_spacing_hz = 1.0 / loop_delay_s
        first_notch_hz = 0.5 / loop_delay_s
        # Notches f_k = (k + 0.5) / tau_loop strictly inside [0, 1000) Hz.
        notches_below_1khz = int(math.floor(1000.0 * loop_delay_s + 0.5))
        icc_cases.append(
            {
                "loop_delay_ms": loop_delay_ms,
                "itu_t_p1150_compliant": bool(loop_delay_ms < icc_p1150_limit_ms),
                "comb_notch_spacing_hz": notch_spacing_hz,
                "first_notch_hz": first_notch_hz,
                "notches_in_0_to_1000_hz": notches_below_1khz,
            }
        )
    comb_peak_to_notch_db = float(20.0 * math.log10((1.0 + icc_alpha) / (1.0 - icc_alpha)))
    closed_loop_gain = 0.5
    howling_margin_db = float(-20.0 * math.log10(closed_loop_gain))

    switching_windows = {
        "window_A_path_change": dual_filter_switch_decision(1.28, 0.08, 0.95),
        "window_B_double_talk_diverged": dual_filter_switch_decision(1.00, 2.50, 0.25),
        "window_C_double_talk_overfit_blocked": dual_filter_switch_decision(1.00, 0.20, 0.30),
    }

    return {
        "stereo_decorrelation": {
            "source": source.tolist(),
            "snapshot_indices": snapshot_indices,
            "true_stereo_path": true_stereo_path.tolist(),
            "lti_allpass": {
                "pole": lti_pole,
                "x2": x2_lti.tolist(),
                "snapshot_matrix": u_lti.tolist(),
                "determinant": float(np.linalg.det(u_lti)),
                "rank": int(np.linalg.matrix_rank(u_lti)),
                "general_lti_rank_upper_bound": 3,
                "null_vector": lti_null_vector.tolist(),
                "null_product": lti_null_residual.tolist(),
            },
            "time_varying_allpass": {
                "pole_trajectory": tv_poles.tolist(),
                "x2": x2_tv.tolist(),
                "snapshot_matrix": u_tv.tolist(),
                "determinant": float(np.linalg.det(u_tv)),
                "rank": int(np.linalg.matrix_rank(u_tv)),
                "echo_snapshots": y_tv.tolist(),
                "recovered_stereo_path": recovered_stereo_path.tolist(),
            },
        },
        "loop_deadlines": {
            "speed_of_sound_mps": c_mps,
            "sample_rate_hz": fs_hz,
            "anc": {
                "primary_distance_m": d_primary_m,
                "secondary_distance_m": d_secondary_m,
                "acoustic_flight_ms": tau_acoust_ms,
                "secondary_path_ms": tau_sec_ms,
                "max_electronic_delay_ms": tau_elec_max_ms,
                "max_electronic_delay_samples": tau_elec_max_samples,
                "pbfdaf_block_samples": pbfdaf_block_samples,
                "pbfdaf_block_ms": pbfdaf_block_ms,
                "block_to_budget_ratio": float(pbfdaf_block_ms / tau_elec_max_ms),
                "broadband_causal_anc_feasible": bool(pbfdaf_block_ms <= tau_elec_max_ms),
            },
            "digital_aec": {
                "system_playback_buffer_samples": aec_buffer_delay_samples,
                "system_playback_buffer_ms": aec_buffer_delay_ms,
                "ram_reference_alignment_feasible": True,
            },
            "icc_itu_t_p1150": {
                "max_reinforcement_delay_ms": icc_p1150_limit_ms,
                "reinforcement_amplitude_ratio": icc_alpha,
                "comb_peak_to_notch_swing_db": comb_peak_to_notch_db,
                "closed_loop_gain": closed_loop_gain,
                "howling_gain_margin_db": howling_margin_db,
                "cases": icc_cases,
            },
        },
        "dual_filter_switching": {
            "copy_ratio": 0.25,
            "reset_ratio": 2.0,
            "coherence_threshold": 0.8,
            "windows": switching_windows,
        },
        "scope": (
            "two-tap deterministic stereo snapshots, free-field flight budget, "
            "two-wave comb model, and window-level dual-filter switching gates"
        ),
    }
