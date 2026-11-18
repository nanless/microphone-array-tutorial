"""E03-19 & E03-20: rigid-sphere baffle scattering bias and MEMS acoustic gasket pole mismatch.

Pure-NumPy deterministic implementations without SciPy or import-time side effects.
Follows the repository forward Fourier kernel exp(-j * 2*pi*f*t) and book azimuth convention.
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


def _spherical_jn_scalar(n: int, x: float) -> float:
    """Evaluate spherical Bessel j_n(x) for integer n >= 0 and x > 0 without upward cancellation."""
    if x <= 10.0:
        term = 1.0
        for m in range(1, n + 1):
            term *= x / float(2 * m + 1)
        acc = term
        half_x2 = 0.5 * x * x
        for k in range(1, 45):
            term *= -half_x2 / float(k * (2 * n + 2 * k + 1))
            acc += term
            if abs(term) <= 1e-17 * abs(acc):
                break
        return acc

    # Downward Miller recurrence for x > 10.0, normalized by exact j_0(x) = sin(x)/x.
    n_start = max(n, int(math.ceil(x))) + 40
    vals = np.zeros(n_start + 2, dtype=np.float64)
    vals[n_start + 1] = 0.0
    vals[n_start] = 1.0
    for m in range(n_start, 0, -1):
        vals[m - 1] = (2.0 * m + 1.0) / x * vals[m] - vals[m + 1]
        if abs(vals[m - 1]) > 1e100:
            vals[: m + 1] *= 1e-100
    j0_exact = math.sin(x) / x
    j1_exact = math.sin(x) / (x * x) - math.cos(x) / x
    scale = (j0_exact / vals[0]) if abs(j0_exact) >= abs(j1_exact) else (j1_exact / vals[1])
    return float(vals[n] * scale)


def _rigid_sphere_mode_weights(ka: float, max_order: int) -> np.ndarray:
    """Compute b_n(ka) = -j / (ka^2 * h_n^(2)'(ka)) for n = 0..max_order under forward FT exp(-j*2*pi*f*t)."""
    jn = np.array([_spherical_jn_scalar(n, ka) for n in range(max_order + 2)], dtype=np.float64)
    yn = np.zeros(max_order + 2, dtype=np.float64)
    yn[0] = -math.cos(ka) / ka
    yn[1] = -math.cos(ka) / (ka * ka) - math.sin(ka) / ka
    overflow_order = max_order + 1
    for m in range(1, max_order + 1):
        nxt = (2.0 * m + 1.0) / ka * yn[m] - yn[m - 1]
        if not math.isfinite(nxt) or abs(nxt) > 1e140:
            overflow_order = m
            break
        yn[m + 1] = nxt

    b_modes = np.zeros(max_order + 1, dtype=np.complex128)
    # n = 0: h_0^(2)'(x) = -h_1^(2)(x) = -(j_1(x) - j*y_1(x))
    dh0 = -(jn[1] - 1j * yn[1])
    b_modes[0] = -1j / (ka * ka * dh0)

    limit_n = min(max_order, overflow_order)
    for n in range(1, limit_n + 1):
        djn = jn[n - 1] - float(n + 1) / ka * jn[n]
        dyn = yn[n - 1] - float(n + 1) / ka * yn[n]
        dhn = djn - 1j * dyn
        denom = ka * ka * dhn
        if math.isfinite(denom.real) and math.isfinite(denom.imag) and abs(denom) > 0.0:
            b_modes[n] = -1j / denom
    return b_modes


def _legendre_table(u: np.ndarray, max_order: int) -> np.ndarray:
    """Evaluate Legendre polynomials P_0(u)..P_max_order(u) with shape (max_order + 1, ...)."""
    out = np.empty((max_order + 1,) + u.shape, dtype=np.float64)
    out[0] = 1.0
    if max_order >= 1:
        out[1] = u
    for n in range(2, max_order + 1):
        out[n] = ((2.0 * n - 1.0) * u * out[n - 1] - float(n - 1) * out[n - 2]) / float(n)
    return out


def rigid_sphere_antipodal_response(
    radius_m: float,
    frequencies_hz: np.ndarray,
    azimuth_rad: np.ndarray,
    *,
    sound_speed: float = SPEED_OF_SOUND,
    max_order: int = 24,
) -> dict[str, np.ndarray]:
    """Exact truncated Helmholtz rigid-sphere surface pressures at (-a,0,0) and (+a,0,0).

    Uses pure-NumPy spherical Bessel/Neumann functions under exp(-j*2*pi*f*t)
    (outgoing Hankel h_n^(2) = j_n - j*y_n) to verify Rayleigh 1.5x dipole limit.
    """
    a_m = finite_real_scalar(radius_m, "radius_m")
    c_m_s = finite_real_scalar(sound_speed, "sound_speed")
    if a_m <= 0.0 or c_m_s <= 0.0:
        raise ValueError("radius_m and sound_speed must be strictly positive")
    if isinstance(max_order, (bool, np.bool_)) or not isinstance(max_order, (int, np.integer)) or max_order < 2:
        raise ValueError("max_order must be an integer >= 2")

    freqs = validate_frequencies(frequencies_hz)
    if np.any(freqs <= 0.0):
        raise ValueError("frequencies_hz must be strictly positive for Helmholtz scattering evaluation")
    az = np.atleast_1d(finite_real_array(azimuth_rad, "azimuth_rad"))
    if az.ndim != 1 or az.size < 1:
        raise ValueError("azimuth_rad must be a scalar or non-empty 1-D real array")

    order = int(max_order)
    u1 = np.sin(az)  # cos(Theta_1) at (+a, 0, 0) for source direction [sin(az), cos(az), 0]
    p_leg = _legendre_table(u1, order)  # shape (order + 1, A)
    n_idx = np.arange(order + 1, dtype=np.float64)
    j_powers = (1j) ** np.arange(order + 1)
    minus_j_powers = (-1j) ** np.arange(order + 1)
    weights_2n1 = 2.0 * n_idx + 1.0

    p0 = np.empty((freqs.size, az.size), dtype=np.complex128)
    p1 = np.empty((freqs.size, az.size), dtype=np.complex128)
    ka_arr = 2.0 * math.pi * freqs * (a_m / c_m_s)

    for idx_f, ka in enumerate(ka_arr):
        b_n = _rigid_sphere_mode_weights(float(ka), order)
        coeff_1 = weights_2n1 * j_powers * b_n
        coeff_0 = weights_2n1 * minus_j_powers * b_n
        p1[idx_f] = coeff_1 @ p_leg
        p0[idx_f] = coeff_0 @ p_leg

    phase_diff_10 = np.angle(p1 / p0)
    tau_10 = -phase_diff_10 / (2.0 * math.pi * freqs[:, np.newaxis])
    tau_free = np.broadcast_to(
        -(2.0 * a_m * np.sin(az)[np.newaxis, :]) / c_m_s,
        (freqs.size, az.size),
    ).copy()

    return {
        "ka": ka_arr,
        "pressure_mic0_left": p0,
        "pressure_mic1_right": p1,
        "phase_diff_10_rad": phase_diff_10,
        "tau_10_sphere_s": tau_10,
        "tau_10_free_s": tau_free,
    }


def _solve_woodworth_overflow_rad() -> float:
    """Solve theta + sin(theta) = 2.0 on (0, pi/2) via Newton's method."""
    theta = 1.106
    for _ in range(15):
        f_val = theta + math.sin(theta) - 2.0
        df_val = 1.0 + math.cos(theta)
        theta -= f_val / df_val
    return theta


def rigid_baffle_aperture_example(
    radius_m: float = 0.05,
    sound_speed: float = SPEED_OF_SOUND,
) -> dict[str, object]:
    """E03-19: Rayleigh 1.5x expansion vs Woodworth diffraction and free-field arcsin DOA bias/overflow."""
    a_m = finite_real_scalar(radius_m, "radius_m")
    c_m_s = finite_real_scalar(sound_speed, "sound_speed")
    if a_m <= 0.0 or c_m_s <= 0.0:
        raise ValueError("radius_m and sound_speed must be strictly positive")

    angles_deg = [20.0, 30.0, 45.0, 60.0]
    rows = []
    for deg in angles_deg:
        rad = math.radians(deg)
        s = math.sin(rad)
        tau_free = -(2.0 * a_m * s) / c_m_s
        tau_lf = 1.5 * tau_free
        tau_hf = -(a_m * (rad + s)) / c_m_s

        arg_lf = -c_m_s * tau_lf / (2.0 * a_m)
        if abs(arg_lf) <= 1.0:
            hat_lf = math.degrees(math.asin(arg_lf))
            bias_lf: float | None = hat_lf - deg
            overflow_lf = False
        else:
            hat_lf = None
            bias_lf = None
            overflow_lf = True

        arg_hf = -c_m_s * tau_hf / (2.0 * a_m)
        if abs(arg_hf) <= 1.0:
            hat_hf = math.degrees(math.asin(arg_hf))
            bias_hf: float | None = hat_hf - deg
            overflow_hf = False
        else:
            hat_hf = None
            bias_hf = None
            overflow_hf = True

        rows.append(
            {
                "true_azimuth_deg": deg,
                "tau_free_us": tau_free * 1e6,
                "tau_rayleigh_lf_us": tau_lf * 1e6,
                "free_field_arcsin_arg_lf": arg_lf,
                "theta_hat_lf_deg": hat_lf,
                "bias_lf_deg": bias_lf,
                "lf_arcsin_overflow": overflow_lf,
                "tau_woodworth_hf_us": tau_hf * 1e6,
                "woodworth_expansion_factor": 0.5 * (1.0 + rad / s),
                "free_field_arcsin_arg_hf": arg_hf,
                "theta_hat_hf_deg": hat_hf,
                "bias_hf_deg": bias_hf,
                "hf_arcsin_overflow": overflow_hf,
            }
        )

    helmholtz_check = rigid_sphere_antipodal_response(
        a_m,
        np.array([20.0, 200.0], dtype=np.float64),
        np.array([math.radians(30.0)], dtype=np.float64),
        sound_speed=c_m_s,
    )
    ratio_20hz = float(helmholtz_check["tau_10_sphere_s"][0, 0] / helmholtz_check["tau_10_free_s"][0, 0])
    ratio_200hz = float(helmholtz_check["tau_10_sphere_s"][1, 0] / helmholtz_check["tau_10_free_s"][1, 0])

    crit_lf_deg = math.degrees(math.asin(2.0 / 3.0))
    crit_hf_deg = math.degrees(_solve_woodworth_overflow_rad())

    return {
        "radius_m": a_m,
        "diameter_m": 2.0 * a_m,
        "sound_speed_m_s": c_m_s,
        "rayleigh_lf_expansion_factor": 1.5,
        "woodworth_expansion_90deg": 0.5 * (1.0 + 0.5 * math.pi),
        "critical_overflow_lf_deg": crit_lf_deg,
        "critical_overflow_hf_deg": crit_hf_deg,
        "helmholtz_exact_ratio_30deg_20hz": ratio_20hz,
        "helmholtz_exact_ratio_30deg_200hz": ratio_200hz,
        "angle_table": rows,
        "limits": (
            "Flush-mounted antipodal sensors on an ideal rigid sphere in free space; "
            "Woodworth creeping-wave delay applies on [-90, 90] deg for ka >> 1."
        ),
    }


def acoustic_vent_highpass(
    frequencies_hz: np.ndarray,
    cutoff_hz: float,
) -> np.ndarray:
    """Evaluate causal first-order acoustic-vent response j*(f/f_c) / (1 + j*(f/f_c))."""
    freqs = validate_frequencies(frequencies_hz)
    f_c = finite_real_scalar(cutoff_hz, "cutoff_hz")
    if f_c <= 0.0:
        raise ValueError("cutoff_hz must be strictly positive")
    ratio = freqs / f_c
    return (1j * ratio) / (1.0 + 1j * ratio)


def mems_gasket_leakage_example(
    sealed_cutoff_hz: float = 25.0,
    leaky_cutoff_hz: float = 80.0,
    spacing_m: float = 0.02,
    sound_speed: float = SPEED_OF_SOUND,
) -> dict[str, object]:
    """E03-20: 1 kHz QA invisibility vs 100 Hz 24.62-deg phase mismatch (11.73x propagation) & DMA inversion."""
    fc1 = finite_real_scalar(sealed_cutoff_hz, "sealed_cutoff_hz")
    fc2 = finite_real_scalar(leaky_cutoff_hz, "leaky_cutoff_hz")
    d_m = finite_real_scalar(spacing_m, "spacing_m")
    c_m_s = finite_real_scalar(sound_speed, "sound_speed")
    if fc1 <= 0.0 or fc2 <= 0.0 or d_m <= 0.0 or c_m_s <= 0.0:
        raise ValueError("cutoff frequencies, spacing_m, and sound_speed must be strictly positive")

    eval_freqs = np.array([100.0, 1000.0], dtype=np.float64)
    h1 = acoustic_vent_highpass(eval_freqs, fc1)
    h2 = acoustic_vent_highpass(eval_freqs, fc2)

    h1_100, h1_1k = complex(h1[0]), complex(h1[1])
    h2_100, h2_1k = complex(h2[0]), complex(h2[1])

    phase_diff_1k_deg = math.degrees(math.atan2((h2_1k / h1_1k).imag, (h2_1k / h1_1k).real))
    level_diff_1k_db = 20.0 * math.log10(abs(h2_1k) / abs(h1_1k))

    phase_1_100_deg = math.degrees(math.atan2(h1_100.imag, h1_100.real))
    phase_2_100_deg = math.degrees(math.atan2(h2_100.imag, h2_100.real))
    delta_phi_vent_rad = math.atan2((h2_100 / h1_100).imag, (h2_100 / h1_100).real)
    delta_phi_vent_deg = math.degrees(delta_phi_vent_rad)
    level_diff_100_db = 20.0 * math.log10(abs(h2_100) / abs(h1_100))

    tau_max = d_m / c_m_s
    delta_phi_prop_rad = 2.0 * math.pi * 100.0 * tau_max
    delta_phi_prop_deg = math.degrees(delta_phi_prop_rad)
    phase_error_to_prop_ratio = delta_phi_vent_rad / delta_phi_prop_rad
    spurious_delay_s = delta_phi_vent_rad / (2.0 * math.pi * 100.0)

    # First-order cardioid DMA at 100 Hz: Y = X_front - exp(-j * delta_phi_prop) * X_back
    # Front wave: X_front = H_1, X_back = H_2 * exp(-j * delta_phi_prop)
    # Back wave:  X_front = H_1 * exp(-j * delta_phi_prop), X_back = H_2
    steer_delay = complex(np.exp(-1j * delta_phi_prop_rad))
    ideal_front_gain = abs(1.0 - steer_delay * steer_delay)
    ideal_back_gain = 0.0
    leaky_front_gain = abs(h1_100 - steer_delay * (h2_100 * steer_delay))
    leaky_back_gain = abs(h1_100 * steer_delay - steer_delay * h2_100)
    leaky_fbr_db = 20.0 * math.log10(leaky_front_gain / leaky_back_gain)

    return {
        "sealed_cutoff_hz": fc1,
        "leaky_cutoff_hz": fc2,
        "spacing_m": d_m,
        "sound_speed_m_s": c_m_s,
        "qa_at_1000hz": {
            "mag_sealed": abs(h1_1k),
            "mag_leaky": abs(h2_1k),
            "relative_level_db": level_diff_1k_db,
            "phase_error_deg": phase_diff_1k_deg,
        },
        "mismatch_at_100hz": {
            "mag_sealed": abs(h1_100),
            "mag_leaky": abs(h2_100),
            "relative_gain_ratio": abs(h2_100) / abs(h1_100),
            "relative_level_db": level_diff_100_db,
            "phase_sealed_deg": phase_1_100_deg,
            "phase_leaky_deg": phase_2_100_deg,
            "phase_error_rad": delta_phi_vent_rad,
            "phase_error_deg": delta_phi_vent_deg,
            "max_propagation_phase_deg": delta_phi_prop_deg,
            "max_propagation_delay_us": tau_max * 1e6,
            "error_to_propagation_phase_ratio": phase_error_to_prop_ratio,
            "spurious_time_shift_us": spurious_delay_s * 1e6,
        },
        "cardioid_dma_100hz": {
            "ideal_front_gain": ideal_front_gain,
            "ideal_back_gain": ideal_back_gain,
            "leaky_front_gain": leaky_front_gain,
            "leaky_back_gain": leaky_back_gain,
            "leaky_front_to_back_db": leaky_fbr_db,
        },
        "limits": (
            "First-order lumped RC acoustic-vent high-pass model; "
            "does not include Helmholtz port resonance at ultrasonic frequencies."
        ),
    }
