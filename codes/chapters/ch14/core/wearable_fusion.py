"""Cross-modal air/in-ear own-voice fusion, hear-through comb response, and wind coherence."""
from __future__ import annotations

import numpy as np

from codes.chapters.ch02.core.conventions import (
    SPEED_OF_SOUND,
    finite_real_scalar,
    hermitian_part,
    validate_frequencies,
)


def _validate_complex(value, name: str) -> np.ndarray:
    original = np.asarray(value)
    if original.dtype.kind not in "iufc":
        raise ValueError(f"{name} must contain numeric values")
    with np.errstate(over="ignore", invalid="ignore"):
        arr = np.asarray(original, dtype=np.complex128)
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be finite")
    return arr


def cross_modal_mvdr_mmse(
    c_occ: complex | float | np.ndarray,
    R_nn: np.ndarray,
    phi_s: float = 1.0,
) -> dict:
    """Compute distortionless MVDR and MMSE weights for cross-modal own-voice fusion.

    Parameters
    ----------
    c_occ : complex scalar ``H_occ`` or 1-D vector ``c`` of shape ``(M,)``
        If scalar, constructs the 2-channel relative transfer vector ``[1, H_occ]^T``
        referencing the outer air microphone ``Y_air``. If 1-D vector, uses ``c_occ``
        directly as the relative transfer vector ``c``.
    R_nn : (M, M) Hermitian positive definite noise covariance matrix.
    phi_s : float, default 1.0
        Clean reference target speech variance ``E[|S_air|^2] > 0``.
    """
    phi_val = finite_real_scalar(phi_s, "phi_s")
    if phi_val <= 0.0:
        raise ValueError("phi_s must be strictly positive")

    r_raw = _validate_complex(R_nn, "R_nn")
    if r_raw.ndim != 2 or r_raw.shape[0] != r_raw.shape[1] or r_raw.shape[0] < 2:
        raise ValueError("R_nn must be a square matrix of dimension at least 2")
    n_chan = r_raw.shape[0]

    scale_r = float(np.max(np.abs(r_raw)))
    if scale_r == 0.0:
        raise ValueError("R_nn must be non-zero and Hermitian")
    r_scaled = r_raw / scale_r
    asym = float(np.linalg.norm(r_scaled - r_scaled.conj().T, ord="fro"))
    norm_r = float(np.linalg.norm(r_scaled, ord="fro"))
    if norm_r == 0.0 or asym > 1e-10 * norm_r:
        raise ValueError("R_nn must be non-zero and Hermitian")
    r_mat = hermitian_part(r_raw)

    c_raw = _validate_complex(c_occ, "c_occ")
    if c_raw.ndim == 0:
        if n_chan != 2:
            raise ValueError("scalar c_occ requires 2x2 R_nn")
        c_vec = np.array([1.0 + 0.0j, c_raw.item()], dtype=np.complex128)
    elif c_raw.ndim == 1 and c_raw.shape == (n_chan,):
        c_vec = c_raw
    else:
        raise ValueError(f"c_occ must be a scalar or 1-D vector of shape ({n_chan},)")

    if float(np.linalg.norm(c_vec)) == 0.0:
        raise ValueError("relative transfer vector c must be non-zero")

    if 1e-100 <= scale_r <= 1e100:
        try:
            chol = np.linalg.cholesky(r_mat)
        except np.linalg.LinAlgError as exc:
            raise ValueError("R_nn must be strictly positive definite") from exc
        y_tmp = np.linalg.solve(chol, c_vec)
        r_inv_c = np.linalg.solve(chol.conj().T, y_tmp)
        q_quad = float(np.real(np.vdot(y_tmp, y_tmp)))
        if not np.isfinite(q_quad) or q_quad <= 0.0:
            raise ValueError("quadratic SNR form c^H R_nn^{-1} c must be finite and positive")
        w_mvdr = r_inv_c / q_quad
        mvdr_noise_var = 1.0 / q_quad
        with np.errstate(over="raise", invalid="raise"):
            try:
                sinr_linear = phi_val * q_quad
            except FloatingPointError as exc:
                raise ValueError("output SINR overflowed float64") from exc
    else:
        r_unit = hermitian_part(r_scaled)
        try:
            chol_u = np.linalg.cholesky(r_unit)
        except np.linalg.LinAlgError as exc:
            raise ValueError("R_nn must be strictly positive definite") from exc
        y_u = np.linalg.solve(chol_u, c_vec)
        r_u_inv_c = np.linalg.solve(chol_u.conj().T, y_u)
        q_u = float(np.real(np.vdot(y_u, y_u)))
        if not np.isfinite(q_u) or q_u <= 0.0:
            raise ValueError("quadratic SNR form c^H R_nn^{-1} c must be finite and positive")
        w_mvdr = r_u_inv_c / q_u
        with np.errstate(over="raise", invalid="raise"):
            try:
                mvdr_noise_var = scale_r / q_u
                sinr_linear = (phi_val / scale_r) * q_u
            except FloatingPointError as exc:
                raise ValueError("cross-modal MVDR variance or SINR overflowed float64") from exc

    if not (
        np.all(np.isfinite(w_mvdr))
        and np.isfinite(mvdr_noise_var)
        and np.isfinite(sinr_linear)
        and sinr_linear > 0.0
    ):
        raise ValueError("cross-modal MVDR solution or SINR is non-finite")
    sinr_db = float(10.0 * np.log10(sinr_linear))

    wiener_gain = sinr_linear / (1.0 + sinr_linear)
    w_mmse = wiener_gain * w_mvdr
    mmse_mse = phi_val / (1.0 + sinr_linear)
    if not (np.isfinite(wiener_gain) and np.all(np.isfinite(w_mmse)) and np.isfinite(mmse_mse)):
        raise ValueError("cross-modal MMSE solution is non-finite")

    air_var = float(np.real(r_mat[0, 0]))
    c0_abs = float(np.abs(c_vec[0]))
    if c0_abs > 0.0:
        air_ratio = (phi_val / air_var) * (c0_abs**2)
        air_sinr_db = float(10.0 * np.log10(air_ratio))
    else:
        air_sinr_db = -np.inf

    in_var = float(np.real(r_mat[-1, -1]))
    c_in_abs = float(np.abs(c_vec[-1]))
    if c_in_abs > 0.0:
        equalized_in_noise_var = (in_var / c_in_abs) / c_in_abs
        in_ratio = (phi_val / in_var) * (c_in_abs**2)
        in_sinr_db = float(10.0 * np.log10(in_ratio))
    else:
        equalized_in_noise_var = np.inf
        in_sinr_db = -np.inf

    return {
        "transfer_vector": c_vec,
        "w_mvdr": w_mvdr,
        "w_mmse": w_mmse,
        "wiener_postfilter_gain": float(wiener_gain),
        "mvdr_noise_variance": float(mvdr_noise_var),
        "mmse_total_mse": float(mmse_mse),
        "output_sinr_linear": float(sinr_linear),
        "output_sinr_db": sinr_db,
        "air_only_sinr_db": air_sinr_db,
        "inear_only_sinr_db": in_sinr_db,
        "equalized_inear_noise_variance": float(equalized_in_noise_var),
    }


def hear_through_response(
    freqs_hz: np.ndarray,
    tau_proc_s: float,
    g_leak: float = 1.0,
    g_ht: float = 1.0,
) -> dict:
    """Compute hear-through superposition transfer function and comb-notch frequencies.

    Models ``H_total(f) = g_leak + g_ht * exp(-1j * 2*pi*f*tau_proc_s)``.
    For ``tau_proc_s > 0`` and ``g_leak * g_ht > 0``, destructive notches occur at
    ``f_m = (2*m + 1) / (2 * tau_proc_s)`` for ``m = 0, 1, 2, ...``.
    """
    freqs = validate_frequencies(freqs_hz)
    tau = finite_real_scalar(tau_proc_s, "tau_proc_s")
    if tau < 0.0:
        raise ValueError("tau_proc_s must be non-negative")
    gleak = finite_real_scalar(g_leak, "g_leak")
    ght = finite_real_scalar(g_ht, "g_ht")

    with np.errstate(over="raise", invalid="raise"):
        try:
            phase = -2.0 * np.pi * freqs * tau
            h_total = gleak + ght * np.exp(1j * phase)
            power_gain = (
                np.float64(gleak) ** 2
                + np.float64(ght) ** 2
                + 2.0 * np.float64(gleak) * np.float64(ght) * np.cos(2.0 * np.pi * freqs * tau)
            )
        except (FloatingPointError, OverflowError) as exc:
            raise ValueError("hear-through response overflowed float64") from exc
    if not (np.all(np.isfinite(h_total)) and np.all(np.isfinite(power_gain))):
        raise ValueError("hear-through response is non-finite")
    power_gain = np.maximum(0.0, power_gain)

    f_max = float(np.max(freqs))
    same_sign = (gleak > 0.0 and ght > 0.0) or (gleak < 0.0 and ght < 0.0)
    opposite_sign = (gleak > 0.0 and ght < 0.0) or (gleak < 0.0 and ght > 0.0)
    if tau > 0.0 and same_sign:
        # (2*m + 1) / (2*tau) <= f_max  =>  m <= floor(f_max * tau - 0.5)
        max_m = int(np.floor(f_max * tau - 0.5))
        if max_m >= 0:
            m_idx = np.arange(max_m + 1, dtype=np.float64)
            notches = (2.0 * m_idx + 1.0) / (2.0 * tau)
        else:
            notches = np.empty(0, dtype=np.float64)
    elif tau > 0.0 and opposite_sign:
        max_m = int(np.floor(f_max * tau))
        notches = np.arange(max_m + 1, dtype=np.float64) / tau
    else:
        notches = np.empty(0, dtype=np.float64)

    return {
        "frequencies_hz": freqs,
        "tau_proc_s": tau,
        "g_leak": gleak,
        "g_ht": ght,
        "transfer_function": h_total,
        "power_gain": power_gain,
        "comb_notch_frequencies_hz": notches,
        "first_notch_hz": float(1.0 / (2.0 * tau)) if (tau > 0.0 and gleak * ght > 0.0) else None,
    }


def corcos_wind_coherence(
    freqs_hz: np.ndarray,
    d_m: float,
    u_conv_mps: float,
    alpha_long: float = 0.14,
    *,
    sound_speed: float = SPEED_OF_SOUND,
) -> dict:
    """Compare Corcos convective wind coherence against 3-D diffuse acoustic coherence."""
    freqs = validate_frequencies(freqs_hz)
    d_val = finite_real_scalar(d_m, "d_m")
    u_val = finite_real_scalar(u_conv_mps, "u_conv_mps")
    alpha_val = finite_real_scalar(alpha_long, "alpha_long")
    c_val = finite_real_scalar(sound_speed, "sound_speed")

    if d_val < 0.0:
        raise ValueError("d_m must be non-negative")
    if u_val <= 0.0 or alpha_val <= 0.0 or c_val <= 0.0:
        raise ValueError("u_conv_mps, alpha_long, and sound_speed must be strictly positive")

    with np.errstate(over="raise", invalid="raise"):
        try:
            kc_d = 2.0 * np.pi * freqs * d_val / u_val
            k_ac_d = 2.0 * np.pi * freqs * d_val / c_val
        except FloatingPointError as exc:
            raise ValueError("wind/acoustic wavenumber product overflowed float64") from exc
    wind_coh = np.exp(-alpha_val * kc_d)
    # np.sinc(x) computes sin(pi*x)/(pi*x), so 2*f*d/c gives sin(2*pi*f*d/c)/(2*pi*f*d/c).
    acoustic_coh = np.sinc(2.0 * freqs * d_val / c_val)
    if not (np.all(np.isfinite(wind_coh)) and np.all(np.isfinite(acoustic_coh))):
        raise ValueError("coherence values must be finite")

    return {
        "frequencies_hz": freqs,
        "d_m": d_val,
        "u_conv_mps": u_val,
        "alpha_long": alpha_val,
        "sound_speed_mps": c_val,
        "convective_kc_d_rad": kc_d,
        "acoustic_k_d_rad": k_ac_d,
        "wind_coherence": wind_coh,
        "acoustic_diffuse_coherence": acoustic_coh,
    }
