"""E02-19 & E02-20: 3D vs 2D diffuse-field coherence and fractional-delay FIR interpolation.

Pure-NumPy deterministic implementations without SciPy or import-time side effects.
"""

from __future__ import annotations

import math

import numpy as np

from codes.chapters.ch02.core.conventions import (
    SPEED_OF_SOUND,
    finite_real_array,
    finite_real_scalar,
    validate_frequencies,
)


FIRST_BESSEL_J0_ROOT = 2.404825557695773


def bessel_j0(x: np.ndarray | float, *, max_terms: int = 80) -> np.ndarray:
    """Evaluate J_0(x) in pure NumPy via convergent power series and Jacobi-Anger/Hankel forms.

    Validates finite real input via ``finite_real_array`` (rejects complex/bool/NaN/Inf).
    Uses exact alternating power series for ``|x| <= 12``, spectrally convergent 128-point
    Jacobi-Anger midpoint quadrature for ``12 < |x| <= 56``, and 7-term Hankel asymptotic
    expansion for ``|x| > 56``, achieving ``< 1e-14`` absolute error for all finite ``float64``.
    """
    if isinstance(max_terms, (bool, np.bool_)) or not isinstance(max_terms, (int, np.integer)) or max_terms < 20:
        raise ValueError("max_terms must be an integer >= 20")
    arr = finite_real_array(x, "x")
    abs_x = np.abs(arr)
    out = np.empty_like(abs_x, dtype=np.float64)

    mask_series = abs_x <= 12.0
    if np.any(mask_series):
        xs = abs_x[mask_series]
        quarter_x2 = 0.25 * xs * xs
        term = np.ones_like(xs, dtype=np.float64)
        acc = term.copy()
        for m in range(1, int(max_terms) + 1):
            term = -term * quarter_x2 / (m * m)
            acc = acc + term
            if float(np.max(np.abs(term))) <= 1e-17 * max(1.0, float(np.max(np.abs(acc)))):
                break
        out[mask_series] = acc

    mask_mid = (abs_x > 12.0) & (abs_x <= 56.0)
    if np.any(mask_mid):
        xm = abs_x[mask_mid]
        # Exact Jacobi-Anger Fourier identity on [0, pi]: 128-pt midpoint alias error < 1e-60 on |x| <= 56.
        theta = (np.arange(128, dtype=np.float64) + 0.5) * (math.pi / 128.0)
        cos_theta = np.cos(theta)
        out[mask_mid] = np.mean(np.cos(xm[..., np.newaxis] * cos_theta), axis=-1)

    mask_asymp = abs_x > 56.0
    if np.any(mask_asymp):
        xa = abs_x[mask_asymp]
        inv8x = 1.0 / (8.0 * xa)
        inv8x2 = inv8x * inv8x
        # Hankel asymptotic series for nu = 0 through c7 (truncation < 1e-17 for |x| > 56).
        c1 = -1.0 * inv8x
        c2 = 4.5 * inv8x2
        c3 = -37.5 * inv8x2 * inv8x
        c4 = (9.0 * 25.0 * 49.0 / 24.0) * inv8x2 * inv8x2
        c5 = -(9.0 * 25.0 * 49.0 * 81.0 / 120.0) * inv8x2 * inv8x2 * inv8x
        c6 = -c5 * (121.0 / 6.0) * inv8x
        c7 = -c6 * (169.0 / 7.0) * inv8x
        p0 = 1.0 - c2 + c4 - c6
        q0 = c1 - c3 + c5 - c7
        chi = xa - 0.25 * math.pi
        out[mask_asymp] = np.sqrt(2.0 / (math.pi * xa)) * (p0 * np.cos(chi) - q0 * np.sin(chi))

    return out


def isotropic_diffuse_coherence_pair(
    frequencies_hz: np.ndarray,
    spacing_m: float,
    *,
    sound_speed: float = SPEED_OF_SOUND,
) -> dict[str, np.ndarray]:
    """Return 3-D spherical sinc(2*f*d/c) and 2-D cylindrical J_0(2*pi*f*d/c) coherence."""
    freqs = validate_frequencies(frequencies_hz)
    d_m = finite_real_scalar(spacing_m, "spacing_m")
    c_m_s = finite_real_scalar(sound_speed, "sound_speed")
    if d_m <= 0.0 or c_m_s <= 0.0:
        raise ValueError("spacing_m and sound_speed must be strictly positive")

    with np.errstate(over="ignore", invalid="ignore"):
        kappa_d = 2.0 * math.pi * freqs * (d_m / c_m_s)
    if not np.all(np.isfinite(kappa_d)):
        raise ValueError("kappa * d exceeded float64 range")

    gamma_3d = np.sinc(2.0 * freqs * (d_m / c_m_s))
    gamma_2d = bessel_j0(kappa_d)
    return {
        "frequencies_hz": freqs,
        "kappa_d_rad": kappa_d,
        "gamma_3d": gamma_3d,
        "gamma_2d": gamma_2d,
    }


def _two_mic_endfire_metrics(kappa_d: float, rho_3d: float, rho_2d: float) -> dict[str, float]:
    """Exact 2-mic endfire superdirective metrics designed for 3D diffuse vs evaluated in 2D."""
    # Steering vector for two-mic endfire: a = [1, exp(+j*kappa_d)]^T (consistent with E02-19).
    a_vec = np.array([1.0 + 0.0j, np.exp(1j * kappa_d)], dtype=np.complex128)
    gamma_3d_mat = np.array([[1.0, rho_3d], [rho_3d, 1.0]], dtype=np.float64)
    gamma_2d_mat = np.array([[1.0, rho_2d], [rho_2d, 1.0]], dtype=np.float64)

    inv_3d_a = np.linalg.solve(gamma_3d_mat, a_vec)
    denom_3d = float(np.vdot(a_vec, inv_3d_a).real)
    w_sd_3d = inv_3d_a / denom_3d

    inv_2d_a = np.linalg.solve(gamma_2d_mat, a_vec)
    denom_2d = float(np.vdot(a_vec, inv_2d_a).real)
    w_sd_2d = inv_2d_a / denom_2d

    p_pred_3d = float(np.vdot(w_sd_3d, gamma_3d_mat @ w_sd_3d).real)
    p_actual_2d = float(np.vdot(w_sd_3d, gamma_2d_mat @ w_sd_3d).real)
    p_opt_2d = float(np.vdot(w_sd_2d, gamma_2d_mat @ w_sd_2d).real)
    wng_lin = 1.0 / float(np.vdot(w_sd_3d, w_sd_3d).real)

    return {
        "kappa_d_rad": kappa_d,
        "gamma_3d": rho_3d,
        "gamma_2d": rho_2d,
        "p_out_3d_predicted": p_pred_3d,
        "di_3d_predicted_db": float(-10.0 * math.log10(p_pred_3d)),
        "p_out_2d_actual": p_actual_2d,
        "di_2d_actual_db": float(-10.0 * math.log10(p_actual_2d)),
        "noise_leakage_2d_vs_3d_db": float(10.0 * math.log10(p_actual_2d / p_pred_3d)),
        "p_out_2d_optimal": p_opt_2d,
        "di_2d_optimal_db": float(-10.0 * math.log10(p_opt_2d)),
        "wng_db": float(10.0 * math.log10(wng_lin)),
    }


def diffuse_field_mismatch_example(
    spacing_m: float = 0.04,
    sound_speed: float = SPEED_OF_SOUND,
) -> dict[str, object]:
    """E02-19: exact zero crossings (4287.5 Hz vs ~3282.0 Hz) and 2-D cabin superdirective leakage."""
    d_m = finite_real_scalar(spacing_m, "spacing_m")
    c_m_s = finite_real_scalar(sound_speed, "sound_speed")
    if d_m <= 0.0 or c_m_s <= 0.0:
        raise ValueError("spacing_m and sound_speed must be strictly positive")

    f_zero_3d = c_m_s / (2.0 * d_m)
    f_zero_2d = c_m_s * FIRST_BESSEL_J0_ROOT / (2.0 * math.pi * d_m)
    gamma_2d_at_3d_zero = float(bessel_j0(math.pi))
    gamma_3d_at_2d_zero = math.sin(FIRST_BESSEL_J0_ROOT) / FIRST_BESSEL_J0_ROOT

    # Asymptotic kappa*d -> 0+ 3D-optimal endfire DMA (alpha = 1/3) placed in 2D cylindrical diffuse field:
    alpha_3d_opt = 1.0 / 3.0
    p_asymp_3d = 0.25
    p_asymp_2d = (alpha_3d_opt * alpha_3d_opt + 0.5) / ((1.0 + alpha_3d_opt) ** 2)  # exact 11/32 = 0.34375

    eval_freqs = np.array([500.0, f_zero_2d, f_zero_3d], dtype=np.float64)
    pair = isotropic_diffuse_coherence_pair(eval_freqs, d_m, sound_speed=c_m_s)
    rows = []
    for f_hz, kd, g3, g2 in zip(
        pair["frequencies_hz"],
        pair["kappa_d_rad"],
        pair["gamma_3d"],
        pair["gamma_2d"],
    ):
        row = {"frequency_hz": float(f_hz)}
        row.update(_two_mic_endfire_metrics(float(kd), float(g3), float(g2)))
        rows.append(row)

    return {
        "spacing_m": d_m,
        "sound_speed_m_s": c_m_s,
        "first_zero_3d_hz": f_zero_3d,
        "first_zero_2d_hz": f_zero_2d,
        "first_zero_2d_relative_drop": (f_zero_3d - f_zero_2d) / f_zero_3d,
        "gamma_2d_at_3d_first_zero": gamma_2d_at_3d_zero,
        "gamma_3d_at_2d_first_zero": gamma_3d_at_2d_zero,
        # Legacy key name retained for compatibility; physically this is the 3D-optimal
        # first-order hypercardioid (DI-maximizing pattern, null at arccos(-1/3) ~= 109.47 deg,
        # parameterized here via (alpha + cos theta) / (1 + alpha) with alpha = 1/3).
        "low_freq_limit_3d_supercardioid_in_2d": {
            "alpha_3d_optimal": alpha_3d_opt,
            "p_out_3d_predicted": p_asymp_3d,
            "di_3d_predicted_db": float(-10.0 * math.log10(p_asymp_3d)),
            "p_out_2d_actual": p_asymp_2d,
            "q_2d_actual": 1.0 / p_asymp_2d,
            "di_2d_actual_db": float(-10.0 * math.log10(p_asymp_2d)),
            "noise_leakage_db": float(10.0 * math.log10(p_asymp_2d / p_asymp_3d)),
        },
        "spot_frequencies": rows,
        "limits": (
            "Narrowband isotropic 3D spherical vs 2D horizontal cylindrical coherence models; "
            "does not model non-isotropic cabin reflections or reverberation time T60."
        ),
    }


def windowed_sinc_fractional_delay(
    delay_samples: float,
    num_taps: int = 8,
    *,
    window: str = "hamming",
    dc_normalize: bool = True,
) -> np.ndarray:
    """Design a causal finite-length windowed-sinc fractional-delay FIR filter."""
    d_samp = finite_real_scalar(delay_samples, "delay_samples")
    if isinstance(num_taps, (bool, np.bool_)) or not isinstance(num_taps, (int, np.integer)) or num_taps < 2:
        raise ValueError("num_taps must be an integer >= 2")
    if not isinstance(dc_normalize, (bool, np.bool_)):
        raise ValueError("dc_normalize must be a bool")
    taps = int(num_taps)
    if not (0.0 <= d_samp <= float(taps - 1)):
        raise ValueError("delay_samples must lie within causal support [0, num_taps - 1]")

    idx = np.arange(taps, dtype=np.float64)
    if window == "hamming":
        win = 0.54 - 0.46 * np.cos(2.0 * math.pi * idx / float(taps - 1))
    elif window == "hann":
        win = 0.5 - 0.5 * np.cos(2.0 * math.pi * idx / float(taps - 1))
    elif window == "rectangular":
        win = np.ones(taps, dtype=np.float64)
    else:
        raise ValueError(f"unsupported window {window!r}; expected 'hamming', 'hann', or 'rectangular'")

    h = win * np.sinc(idx - d_samp)
    if dc_normalize:
        dc_sum = float(np.sum(h))
        if abs(dc_sum) <= 1e-12:
            raise ValueError("windowed-sinc filter has near-zero DC gain and cannot be DC-normalized")
        h = h / dc_sum
    return h


def fractional_delay_notch_example(
    sample_rate_hz: float = 16000.0,
) -> dict[str, object]:
    """E02-20: 4-point linear half-sample calculation, mixed-channel DSB/null loss, and 8/16-tap FIR."""
    fs = finite_real_scalar(sample_rate_hz, "sample_rate_hz")
    if fs <= 0.0:
        raise ValueError("sample_rate_hz must be strictly positive")

    # 1. Exact 4-point DFT & steady-state circular convolution check of h_lin = [0.5, 0.5, 0, 0]
    h_lin_4 = np.array([0.5, 0.5, 0.0, 0.0], dtype=np.float64)
    dft_lin_4 = np.fft.fft(h_lin_4)
    x_dc = np.array([1.0, 1.0, 1.0, 1.0], dtype=np.float64)
    x_quarter = np.array([1.0, 0.0, -1.0, 0.0], dtype=np.float64)
    x_nyq = np.array([1.0, -1.0, 1.0, -1.0], dtype=np.float64)
    y_dc = 0.5 * x_dc + 0.5 * np.roll(x_dc, 1)
    y_quarter = 0.5 * x_quarter + 0.5 * np.roll(x_quarter, 1)
    y_nyq = 0.5 * x_nyq + 0.5 * np.roll(x_nyq, 1)

    # 2. Mixed integer + half-sample two-channel alignment across 1, 4, 6, 8 kHz (scaled with fs)
    frac_ratios = np.array([1.0 / 16.0, 0.25, 0.375, 0.5], dtype=np.float64)
    eval_freqs = frac_ratios * fs
    omega_d = 2.0 * math.pi * frac_ratios

    h_wsinc_8 = windowed_sinc_fractional_delay(3.5, num_taps=8, window="hamming", dc_normalize=True)
    h_wsinc_16 = windowed_sinc_fractional_delay(7.5, num_taps=16, window="hamming", dc_normalize=True)

    rows = []
    for f_hz, om in zip(eval_freqs, omega_d):
        # Linear half-sample aligned relative to exact half-sample phase exp(-j * om * 0.5):
        mag_lin = abs(math.cos(0.5 * float(om)))
        dsb_lin = 0.5 * (1.0 + mag_lin)
        null_lin = abs(1.0 - mag_lin)

        # Causal 8-tap (center delay 3.5 vs integer reference 3 -> relative phase compensated by exp(+j*3.5*om))
        n8 = np.arange(8, dtype=np.float64)
        resp_8 = complex(np.sum(h_wsinc_8 * np.exp(-1j * float(om) * n8)) * np.exp(1j * 3.5 * float(om)))
        dsb_8 = 0.5 * abs(1.0 + resp_8)
        null_8 = abs(1.0 - resp_8)

        n16 = np.arange(16, dtype=np.float64)
        resp_16 = complex(np.sum(h_wsinc_16 * np.exp(-1j * float(om) * n16)) * np.exp(1j * 7.5 * float(om)))
        dsb_16 = 0.5 * abs(1.0 + resp_16)
        null_16 = abs(1.0 - resp_16)

        rows.append(
            {
                "frequency_hz": float(f_hz),
                "omega_d_rad": float(om),
                "linear_mag": mag_lin,
                "linear_mag_db": float(20.0 * math.log10(mag_lin)) if mag_lin > 1e-15 else None,
                "linear_dsb_gain": dsb_lin,
                "linear_dsb_gain_db": float(20.0 * math.log10(dsb_lin)),
                "linear_null_residual": null_lin,
                "linear_null_residual_db": float(20.0 * math.log10(null_lin)) if null_lin > 1e-15 else None,
                "wsinc8_mag": abs(resp_8),
                "wsinc8_dsb_gain": dsb_8,
                "wsinc8_dsb_gain_db": float(20.0 * math.log10(dsb_8)),
                "wsinc8_null_residual": null_8,
                "wsinc8_null_residual_db": float(20.0 * math.log10(null_8)) if null_8 > 1e-15 else None,
                "wsinc16_mag": abs(resp_16),
                "wsinc16_dsb_gain": dsb_16,
                "wsinc16_null_residual": null_16,
                "wsinc16_null_residual_db": float(20.0 * math.log10(null_16)) if null_16 > 1e-15 else None,
            }
        )

    return {
        "sample_rate_hz": fs,
        "four_point_linear_half_sample": {
            "h_lin_4": h_lin_4.tolist(),
            "dft_real": dft_lin_4.real.tolist(),
            "dft_imag": dft_lin_4.imag.tolist(),
            "dft_magnitude": np.abs(dft_lin_4).tolist(),
            "y_dc": y_dc.tolist(),
            "y_quarter_fs": y_quarter.tolist(),
            "y_nyquist": y_nyq.tolist(),
        },
        "wsinc8_coefficients": h_wsinc_8.tolist(),
        "wsinc16_coefficients": h_wsinc_16.tolist(),
        "frequency_comparison": rows,
        "limits": (
            "Even-length symmetric Type-II FIR filters still have an exact structural zero at omega_d = pi "
            "for half-sample delay eta = 0.5; practical systems maintain an anti-aliasing guard band below Nyquist."
        ),
    }
