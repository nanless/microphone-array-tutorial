"""Pure-NumPy teaching kernels for online RLS-WPE silence gating and AEC+WPE cascades.

Covers Chapter 07 exercises E07-25 (online RLS-WPE inverse-correlation exponential
inflation over 1000 silence frames, two-frame onset pop, and convex diagonal
pullback gating) and E07-26 (Linear AEC -> WPE exact identifiability vs
pre-WPE NLP non-commutativity counterexample vs WPE -> AEC reversal, plus
frequency-dependent prediction tap budget K(f)). No SciPy imports; float64/complex128.
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


def _checked_hpd_matrix(matrix, name: str = "matrix") -> np.ndarray:
    raw = np.asarray(matrix)
    if raw.dtype.kind not in "iufc":
        raise ValueError(f"{name} must be numeric")
    arr = np.asarray(raw, dtype=np.complex128 if raw.dtype.kind == "c" else np.float64)
    if arr.ndim != 2 or arr.shape[0] == 0 or arr.shape[0] != arr.shape[1]:
        raise ValueError(f"{name} must be a non-empty square matrix")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be finite")
    if not np.allclose(arr, arr.conj().T, rtol=1e-12, atol=1e-15):
        raise ValueError(f"{name} must be Hermitian")
    try:
        np.linalg.cholesky(arr)
    except np.linalg.LinAlgError as error:
        raise ValueError(f"{name} must be positive definite") from error
    return arr


def convex_diagonal_pullback_step(
    inverse_correlation: np.ndarray,
    *,
    beta: float,
    delta: float,
) -> np.ndarray:
    """Apply gated convex pullback ``P_t = (1 - beta) P_{t-1} + beta * delta^{-1} I_D``."""
    p_prev = _checked_hpd_matrix(inverse_correlation, "inverse_correlation")
    beta_val = _finite_real_scalar(beta, "beta")
    delta_val = _finite_real_scalar(delta, "delta")
    if not 0.0 < beta_val < 1.0:
        raise ValueError("beta must lie strictly in (0, 1)")
    if delta_val <= 0.0:
        raise ValueError("delta must be strictly positive")
    dim = p_prev.shape[0]
    target_diag = (1.0 / delta_val) * np.eye(dim, dtype=p_prev.dtype)
    p_next = (1.0 - beta_val) * p_prev + beta_val * target_diag
    if not np.all(np.isfinite(p_next)):
        raise ValueError("pullback step overflowed float64")
    return p_next


def rls_wpe_single_bin_step(
    weights: np.ndarray,
    inverse_correlation: np.ndarray,
    history: np.ndarray,
    observation: float | complex,
    *,
    alpha: float,
    power: float,
) -> dict:
    """Execute one causal single-bin RLS-WPE step (Eq. 7-18 and SMW recursion)."""
    p_prev = _checked_hpd_matrix(inverse_correlation, "inverse_correlation")
    dim = p_prev.shape[0]
    w_raw = np.asarray(weights)
    q_raw = np.asarray(history)
    if w_raw.dtype.kind not in "iufc" or q_raw.dtype.kind not in "iufc":
        raise ValueError("weights and history must be numeric")
    if isinstance(observation, (bool, np.bool_)) or np.asarray(observation).dtype.kind not in "iufc":
        raise ValueError("observation must be a numeric scalar")
    obs_arr = np.asarray(observation)
    if obs_arr.ndim != 0 or not np.all(np.isfinite(obs_arr)):
        raise ValueError("observation must be a finite scalar")
    alpha_val = _finite_real_scalar(alpha, "alpha")
    power_val = _finite_real_scalar(power, "power")
    if not 0.0 < alpha_val < 1.0:
        raise ValueError("alpha must lie strictly in (0, 1)")
    if power_val <= 0.0:
        raise ValueError("power must be strictly positive")

    is_complex = any(arr.dtype.kind == "c" for arr in (p_prev, w_raw, q_raw, obs_arr))
    dtype = np.complex128 if is_complex else np.float64
    g_prev = np.asarray(w_raw, dtype=dtype)
    q_vec = np.asarray(q_raw, dtype=dtype)
    p_mat = np.asarray(p_prev, dtype=dtype)
    x_obs = complex(obs_arr) if is_complex else float(obs_arr)
    if g_prev.shape != (dim,) or q_vec.shape != (dim,):
        raise ValueError("weights and history must have shape (D,) matching inverse_correlation")
    if not np.all(np.isfinite(g_prev)) or not np.all(np.isfinite(q_vec)):
        raise ValueError("weights and history must be finite")

    with np.errstate(over="raise", invalid="raise", divide="raise"):
        try:
            prior_output = x_obs - (np.vdot(g_prev, q_vec) if is_complex else float(g_prev @ q_vec))
            pq = p_mat @ q_vec
            quad = float(np.vdot(q_vec, pq).real) if is_complex else float(q_vec @ pq)
            denom = alpha_val * power_val + quad
            gain = pq / denom
            qh_p = q_vec.conj() @ p_mat if is_complex else q_vec @ p_mat
            p_next = (p_mat - np.outer(gain, qh_p)) / alpha_val
            p_next = 0.5 * (p_next + p_next.conj().T)
            err_conj = prior_output.conjugate() if is_complex else prior_output
            g_next = g_prev + gain * err_conj
        except FloatingPointError as error:
            raise ValueError("RLS-WPE update exceeded float64 range") from error

    if not (np.all(np.isfinite(p_next)) and np.all(np.isfinite(g_next)) and np.isfinite(prior_output)):
        raise ValueError("RLS-WPE update produced non-finite values")
    return {
        "prior_output": prior_output,
        "gain": gain,
        "posterior_weights": g_next,
        "posterior_inverse_correlation": p_next,
    }


def rls_wpe_silence_pullback() -> dict:
    """E07-25: 1000-frame silence inflation & first-syllable pop vs convex diagonal pullback."""
    dim = 2
    alpha = 0.99
    silence_frames = 1000
    hop_ms = 8.0
    p0 = np.eye(dim, dtype=np.float64)
    g0 = np.array([0.5, 0.0], dtype=np.float64)

    # 1. Unguarded silence inflation P_{1000} = alpha^{-1000} I_2 and 2-frame onset pop.
    inflation_scale = float(alpha ** (-silence_frames))
    p_unguarded = inflation_scale * p0
    q_1001 = np.array([0.01, 0.0], dtype=np.float64)
    x_1001 = 1.0
    q_1002 = np.array([1.0, 0.0], dtype=np.float64)
    x_1002 = 1.0

    u_step1 = rls_wpe_single_bin_step(
        g0, p_unguarded, q_1001, x_1001, alpha=alpha, power=1.0
    )
    u_step2 = rls_wpe_single_bin_step(
        u_step1["posterior_weights"],
        u_step1["posterior_inverse_correlation"],
        q_1002,
        x_1002,
        alpha=alpha,
        power=1.0,
    )
    unguarded_pop_power = float(u_step2["prior_output"] ** 2)
    unguarded_pop_db = float(10.0 * math.log10(unguarded_pop_power))

    # 2. Anisotropic 2-step exact rational hand check of Eq. (7-36).
    p_aniso_0 = np.diag([4.0, 2.0])
    p_aniso_1 = convex_diagonal_pullback_step(p_aniso_0, beta=0.25, delta=0.5)
    p_aniso_2 = convex_diagonal_pullback_step(p_aniso_1, beta=0.25, delta=0.5)

    # 3. Guarded 1000-frame silence trajectory with beta=0.05, delta=1.0 on P_0 = I_2.
    p_guarded = p0.copy()
    for _ in range(silence_frames):
        p_guarded = convex_diagonal_pullback_step(p_guarded, beta=0.05, delta=1.0)
    g_step1 = rls_wpe_single_bin_step(
        g0, p_guarded, q_1001, x_1001, alpha=alpha, power=1.0
    )
    g_step2 = rls_wpe_single_bin_step(
        g_step1["posterior_weights"],
        g_step1["posterior_inverse_correlation"],
        q_1002,
        x_1002,
        alpha=alpha,
        power=1.0,
    )
    guarded_step2_power = float(g_step2["prior_output"] ** 2)

    return {
        "dimension": dim,
        "alpha": alpha,
        "silence_frames": silence_frames,
        "hop_ms": hop_ms,
        "silence_duration_s": float(silence_frames * hop_ms * 1e-3),
        "initial_weights": g0.tolist(),
        "initial_inverse_correlation": p0.tolist(),
        "unguarded": {
            "inflation_factor": inflation_scale,
            "post_silence_inverse_correlation": p_unguarded.tolist(),
            "post_silence_trace": float(np.trace(p_unguarded)),
            "frame_1001_prior_output": float(u_step1["prior_output"]),
            "frame_1001_gain": u_step1["gain"].tolist(),
            "frame_1001_posterior_weights": u_step1["posterior_weights"].tolist(),
            "frame_1002_prior_output": float(u_step2["prior_output"]),
            "frame_1002_power": unguarded_pop_power,
            "frame_1002_overshoot_db": unguarded_pop_db,
        },
        "anisotropic_pullback_check": {
            "beta": 0.25,
            "delta": 0.5,
            "p0": p_aniso_0.tolist(),
            "p1": p_aniso_1.tolist(),
            "p2": p_aniso_2.tolist(),
            "trace_p0": float(np.trace(p_aniso_0)),
            "trace_p1": float(np.trace(p_aniso_1)),
            "trace_p2": float(np.trace(p_aniso_2)),
            "uniform_trace_upper_bound": float(max(np.trace(p_aniso_0), dim / 0.5)),
        },
        "guarded_pullback": {
            "beta": 0.05,
            "delta": 1.0,
            "post_silence_inverse_correlation": p_guarded.tolist(),
            "post_silence_trace": float(np.trace(p_guarded)),
            "frame_1001_prior_output": float(g_step1["prior_output"]),
            "frame_1001_gain": g_step1["gain"].tolist(),
            "frame_1001_posterior_weights": g_step1["posterior_weights"].tolist(),
            "frame_1002_prior_output": float(g_step2["prior_output"]),
            "frame_1002_power": guarded_step2_power,
        },
        "scope": (
            "single-bin real D=2 RLS-WPE state control with unit onset power; "
            "history ring buffer advances through silence while P_t gating is compared"
        ),
    }


def frequency_dependent_wpe_budget(
    channels: int,
    band_bin_counts: list[int] | tuple[int, ...],
    band_taps: list[int] | tuple[int, ...],
    *,
    complex_bytes: int = 16,
) -> dict:
    """Compute exact correlation memory and O(D^3) Cholesky proxy for subband K(f) (Eq. 7-38)."""
    m_ch = _integer(channels, "channels", minimum=1)
    c_bytes = _integer(complex_bytes, "complex_bytes", minimum=1)
    if (
        isinstance(band_bin_counts, (str, bytes))
        or isinstance(band_taps, (str, bytes))
        or not hasattr(band_bin_counts, "__iter__")
        or not hasattr(band_taps, "__iter__")
    ):
        raise ValueError("band_bin_counts and band_taps must be iterables of positive integers")
    counts = [_integer(v, "band_bin_count", minimum=1) for v in band_bin_counts]
    taps = [_integer(v, "band_tap", minimum=1) for v in band_taps]
    if not counts or len(counts) != len(taps):
        raise ValueError("band_bin_counts and band_taps must have equal non-empty length")

    total_bins = sum(counts)
    total_corr_bytes = 0
    total_pred_bytes = 0
    total_chol_proxy = 0
    bands = []
    for f_count, k_tap in zip(counts, taps):
        dim = m_ch * k_tap
        corr_b = c_bytes * f_count * (dim ** 2)
        pred_b = c_bytes * f_count * dim * m_ch
        chol_p = f_count * (dim ** 3)
        total_corr_bytes += corr_b
        total_pred_bytes += pred_b
        total_chol_proxy += chol_p
        bands.append(
            {
                "bins": f_count,
                "taps": k_tap,
                "history_dimension": dim,
                "correlation_bytes": corr_b,
                "cholesky_cubic_proxy": chol_p,
            }
        )
    return {
        "channels": m_ch,
        "total_bins": total_bins,
        "complex_bytes": c_bytes,
        "bands": bands,
        "correlation_bytes": total_corr_bytes,
        "correlation_MiB": float(total_corr_bytes / (2 ** 20)),
        "predictor_bytes": total_pred_bytes,
        "cholesky_cubic_proxy": total_chol_proxy,
    }


def _fit_two_tap_tail_wpe(signal: np.ndarray, *, delay: int = 2, tail_indices: tuple[int, ...] = (3, 4, 5)) -> dict:
    sig = _finite_real_array(signal, "signal")
    q_rows = np.array(
        [[sig[t - delay], sig[t - delay - 1]] for t in tail_indices],
        dtype=np.float64,
    )
    target = np.array([sig[t] for t in tail_indices], dtype=np.float64)
    r_mat = q_rows.T @ q_rows
    r_vec = q_rows.T @ target
    det_r = float(np.linalg.det(r_mat))
    g_hat = np.linalg.solve(r_mat, r_vec)
    tail_residual = target - q_rows @ g_hat
    return {
        "input_signal": sig.tolist(),
        "history_rows": q_rows.tolist(),
        "tail_targets": target.tolist(),
        "normal_matrix": r_mat.tolist(),
        "determinant": det_r,
        "cross_vector": r_vec.tolist(),
        "estimated_weights": g_hat.tolist(),
        "tail_residual": tail_residual.tolist(),
        "tail_residual_sum_of_squares": float(tail_residual @ tail_residual),
    }


def cascaded_aec_wpe_ordering() -> dict:
    """E07-26: Linear AEC -> WPE vs pre-WPE NLP mask vs WPE -> AEC, plus K(f) budget."""
    clean_early = np.array([1.0, 1.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)
    true_g = np.array([0.5, 0.25], dtype=np.float64)
    z_rev = np.zeros(6, dtype=np.float64)
    for t in range(6):
        z_rev[t] = clean_early[t]
        if t >= 2:
            z_rev[t] += true_g[0] * z_rev[t - 2]
        if t >= 3:
            z_rev[t] += true_g[1] * z_rev[t - 3]

    far_end = np.array([0.0, 0.0, 2.0, 0.0, 0.0, 0.0], dtype=np.float64)
    echo_path = np.array([1.0, 0.5], dtype=np.float64)
    echo = np.convolve(far_end, echo_path)[:6]
    mic = z_rev + echo

    # Case A: Linear AEC -> WPE (exact linear echo cancellation recovers Z_t).
    linear_aec_out = mic - echo
    case_linear = _fit_two_tap_tail_wpe(linear_aec_out)

    # Case B: Linear AEC -> NLP -> WPE (time-varying mask zeroes frame t=2).
    nlp_mask = np.array([1.0, 1.0, 0.0, 1.0, 1.0, 1.0], dtype=np.float64)
    nlp_corrupted = nlp_mask * linear_aec_out
    case_nlp = _fit_two_tap_tail_wpe(nlp_corrupted)
    # Non-commutativity defect at t=4,5: g_t * (G * Z)_t - G * (g * Z)_t.
    commutator_defect = {}
    for t in (4, 5):
        lhs = nlp_mask[t] * (true_g[0] * z_rev[t - 2] + true_g[1] * z_rev[t - 3])
        rhs = true_g[0] * nlp_corrupted[t - 2] + true_g[1] * nlp_corrupted[t - 3]
        commutator_defect[str(t)] = float(lhs - rhs)

    # Case C: WPE -> Linear AEC (unremoved echo burst hijacks WPE normal equations).
    case_reversed = _fit_two_tap_tail_wpe(mic)

    # Part 2: Frequency-dependent taps K(f) resource budget (M=8, F=257).
    uniform_budget = frequency_dependent_wpe_budget(8, [257], [12])
    banded_budget = frequency_dependent_wpe_budget(8, [65, 64, 128], [12, 8, 4])
    mem_ratio = float(banded_budget["correlation_bytes"] / uniform_budget["correlation_bytes"])
    chol_ratio = float(banded_budget["cholesky_cubic_proxy"] / uniform_budget["cholesky_cubic_proxy"])

    return {
        "delay": 2,
        "taps": 2,
        "tail_fit_frames": [3, 4, 5],
        "clean_early_speech": clean_early.tolist(),
        "true_wpe_weights": true_g.tolist(),
        "reverberant_near_end": z_rev.tolist(),
        "far_end_reference": far_end.tolist(),
        "true_echo_path": echo_path.tolist(),
        "echo_signal": echo.tolist(),
        "microphone_signal": mic.tolist(),
        "nlp_mask": nlp_mask.tolist(),
        "non_commutativity_defect_frames_4_5": commutator_defect,
        "architectures": {
            "linear_aec_then_wpe": case_linear,
            "linear_aec_then_nlp_then_wpe": case_nlp,
            "wpe_then_linear_aec": case_reversed,
        },
        "frequency_dependent_taps_budget": {
            "uniform_K12": uniform_budget,
            "three_band_K12_8_4": banded_budget,
            "correlation_memory_ratio": mem_ratio,
            "correlation_memory_reduction_percent": float(100.0 * (1.0 - mem_ratio)),
            "cholesky_proxy_ratio": chol_ratio,
            "cholesky_proxy_reduction_percent": float(100.0 * (1.0 - chol_ratio)),
        },
        "scope": (
            "six-frame deterministic single-bin AR reverberation + 2-tap echo control, "
            "plus exact array/cubic-order resource accounting without timing claims"
        ),
    }
