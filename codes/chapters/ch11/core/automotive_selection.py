"""Deterministic smart-cockpit multi-zone, ICC, and conferencing selection helpers.

Reuses Chapter 3 near-field steering, Chapter 5 DSB/LCMV solvers, and Chapter 11
interval verdicts under explicitly declared single-frequency / delay models.
No cabin recording, neural speech model, or certification lab test is executed.
"""

from __future__ import annotations

from collections.abc import Sequence
from fractions import Fraction
import math
from numbers import Integral

import numpy as np

from codes.chapters.ch02.core.conventions import SPEED_OF_SOUND, finite_real_scalar
from codes.chapters.ch03.core.geometry import near_field_steering
from codes.chapters.ch05.core.beamforming import dsb_weights, lcmv_weights
from codes.chapters.ch11.core.selection import upper_limit_verdict


def _complex_pairs(values: np.ndarray | complex) -> list[object]:
    arr = np.asarray(values, dtype=complex)
    return np.stack((arr.real, arr.imag), axis=-1).tolist()


def _exact_decimal_fraction(value: object, name: str) -> Fraction:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must not be a boolean")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, Integral):
        return Fraction(int(value))
    if isinstance(value, (float, np.floating)):
        if isinstance(value, np.floating) and value.dtype.itemsize > 8:
            raise ValueError(f"{name} extended-precision floats are not supported")
        scalar = float(value)
        if not math.isfinite(scalar):
            raise ValueError(f"{name} must be finite")
        return Fraction(str(scalar))
    raise ValueError(f"{name} must be a finite int, float, or Fraction")


def bus_clock_drift_metrics(
    sro_ppm: float,
    duration_s: float,
    eval_freq_hz: float,
    max_phase_error_deg: float,
) -> dict[str, float | str]:
    """Evaluate cross-node sampling-rate offset phase drift against a hard bound."""
    ppm = finite_real_scalar(sro_ppm, "sro_ppm")
    t_s = finite_real_scalar(duration_s, "duration_s")
    f_hz = finite_real_scalar(eval_freq_hz, "eval_freq_hz")
    phi_lim_deg = finite_real_scalar(max_phase_error_deg, "max_phase_error_deg")
    if t_s <= 0.0 or f_hz <= 0.0 or phi_lim_deg <= 0.0:
        raise ValueError("duration_s, eval_freq_hz, and max_phase_error_deg must be positive")
    epsilon = ppm * 1e-6
    if 1.0 + epsilon <= 0.0:
        raise ValueError("1 + sro_ppm * 1e-6 must be strictly positive")

    exact_drift_s = abs(t_s * epsilon / (1.0 + epsilon))
    max_drift_s = phi_lim_deg / (360.0 * f_hz)
    phase_error_deg = 360.0 * f_hz * exact_drift_s
    verdict = upper_limit_verdict([[phase_error_deg, phase_error_deg]], [phi_lim_deg])
    return {
        "sro_ppm": ppm,
        "duration_s": t_s,
        "eval_freq_hz": f_hz,
        "exact_drift_s": exact_drift_s,
        "exact_drift_us": exact_drift_s * 1e6,
        "max_allowed_drift_us": max_drift_s * 1e6,
        "phase_error_deg": phase_error_deg,
        "max_phase_error_deg": phi_lim_deg,
        "verdict": verdict,
    }


def _physical_point_transfer(
    positions: np.ndarray,
    source: np.ndarray,
    frequency_hz: float,
    nominal_distance_m: float,
    sound_speed: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Factor physical point-source transfer via ch03 near_field_steering."""
    steering_rel0 = near_field_steering(
        positions, source, [frequency_hz], reference=0, sound_speed=sound_speed
    )[0]
    distances = np.linalg.norm(source[None, :] - positions, axis=1)
    ref_prop = (nominal_distance_m / distances[0]) * np.exp(
        -2.0j * np.pi * frequency_hz * distances[0] / sound_speed
    )
    transfer = ref_prop * steering_rel0
    return transfer, steering_rel0, distances


def multizone_cabin_selection_case() -> dict[str, object]:
    """Return the E11-28 smart-cockpit bus sync and 2-zone isolation selection case."""
    frequency_hz = 1000.0
    sound_speed = SPEED_OF_SOUND
    nominal_mouth_m = 0.20
    sensor_noise_power = 1e-3
    scattering_pert = 0.05
    calibration_retention = 0.95

    lcr_min_db = 25.0
    wng_min_db = -3.0
    dt_loss_max_db = 1.0
    upper_limits = [-lcr_min_db, -wng_min_db, dt_loss_max_db]

    s_driver = np.array([-0.35, 0.0])
    s_passenger = np.array([0.35, 0.0])

    bus_clocks = {
        "a2b_sync": bus_clock_drift_metrics(0.0, 60.0, 4000.0, 5.0),
        "bt_async_uncorrected": bus_clock_drift_metrics(50.0, 60.0, 4000.0, 5.0),
    }
    bus_verdicts = {k: str(v["verdict"]) for k, v in bus_clocks.items()}

    # Option 1: Central ceiling 2-mic array
    mics_ceiling = np.array([[-0.03, 0.40], [0.03, 0.40]])
    h_d_ceil, a_d_ceil, r_d_ceil = _physical_point_transfer(
        mics_ceiling, s_driver, frequency_hz, nominal_mouth_m, sound_speed
    )
    h_p_ceil, a_p_ceil, r_p_ceil = _physical_point_transfer(
        mics_ceiling, s_passenger, frequency_hz, nominal_mouth_m, sound_speed
    )
    w_ceil = dsb_weights(a_d_ceil)
    p_d_ceil = float(abs(np.vdot(w_ceil, h_d_ceil)) ** 2)
    p_p_ceil = float(abs(np.vdot(w_ceil, h_p_ceil)) ** 2)
    lcr_ceil = 10.0 * math.log10(p_d_ceil / p_p_ceil)
    norm2_ceil = float(np.vdot(w_ceil, w_ceil).real)
    wng_ceil = 10.0 * math.log10(float(abs(np.vdot(w_ceil, a_d_ceil)) ** 2) / norm2_ceil)
    verdict_ceil = upper_limit_verdict(
        [[-lcr_ceil, -lcr_ceil], [-wng_ceil, -wng_ceil], [0.0, 0.0]], upper_limits
    )

    # Option 2: Driver headrest 2-mic array (4 cm spacing)
    mics_hr1 = np.array([[-0.35, -0.20], [-0.35, -0.24]])
    h_d_hr1, a_d_hr1, r_d_hr1 = _physical_point_transfer(
        mics_hr1, s_driver, frequency_hz, nominal_mouth_m, sound_speed
    )
    h_p_hr1, a_p_hr1, r_p_hr1 = _physical_point_transfer(
        mics_hr1, s_passenger, frequency_hz, nominal_mouth_m, sound_speed
    )

    # 2a: Headrest DSB
    w_hr_dsb = dsb_weights(a_d_hr1)
    y_d_pri = complex(np.vdot(w_hr_dsb, h_d_hr1))
    y_p_pri = complex(np.vdot(w_hr_dsb, h_p_hr1))
    p_d_hr_dsb = float(abs(y_d_pri) ** 2)
    p_p_hr_dsb = float(abs(y_p_pri) ** 2)
    lcr_hr_dsb = 10.0 * math.log10(p_d_hr_dsb / p_p_hr_dsb)
    norm2_hr_dsb = float(np.vdot(w_hr_dsb, w_hr_dsb).real)
    wng_hr_dsb = 10.0 * math.log10(float(abs(np.vdot(w_hr_dsb, a_d_hr1)) ** 2) / norm2_hr_dsb)
    verdict_hr_dsb = upper_limit_verdict(
        [[-lcr_hr_dsb, -lcr_hr_dsb], [-wng_hr_dsb, -wng_hr_dsb], [0.0, 0.0]], upper_limits
    )

    # 2b: Headrest 2-mic LCMV nulling nominal Passenger steering vector
    constraints_hr = np.column_stack((a_d_hr1, a_p_hr1))
    w_hr_lcmv = lcmv_weights(np.eye(2), constraints_hr, np.array([1.0, 0.0]))
    norm2_hr_lcmv = float(np.vdot(w_hr_lcmv, w_hr_lcmv).real)
    wng_hr_lcmv = 10.0 * math.log10(float(abs(np.vdot(w_hr_lcmv, a_d_hr1)) ** 2) / norm2_hr_lcmv)
    h_p_hr1_pert = np.array([h_p_hr1[0], h_p_hr1[1] * (1.0 + scattering_pert)], dtype=complex)
    p_d_hr_lcmv = float(abs(np.vdot(w_hr_lcmv, h_d_hr1)) ** 2)
    p_p_hr_lcmv_ideal = float(abs(np.vdot(w_hr_lcmv, h_p_hr1)) ** 2)
    p_p_hr_lcmv_pert = float(abs(np.vdot(w_hr_lcmv, h_p_hr1_pert)) ** 2)
    lcr_hr_lcmv = 10.0 * math.log10(p_d_hr_lcmv / p_p_hr_lcmv_pert)
    verdict_hr_lcmv = upper_limit_verdict(
        [[-lcr_hr_lcmv, -lcr_hr_lcmv], [-wng_hr_lcmv, -wng_hr_lcmv], [0.0, 0.0]], upper_limits
    )

    # Passenger headrest reference mic 0 at (+0.35, -0.20) on A2B synchronous bus
    mics_hr2 = np.array([[0.35, -0.20], [0.35, -0.24]])
    h_d_hr2, _, r_d_hr2 = _physical_point_transfer(
        mics_hr2, s_driver, frequency_hz, nominal_mouth_m, sound_speed
    )
    h_p_hr2, _, r_p_hr2 = _physical_point_transfer(
        mics_hr2, s_passenger, frequency_hz, nominal_mouth_m, sound_speed
    )
    h_d_ref = complex(h_d_hr2[0])
    h_p_ref = complex(h_p_hr2[0])

    g_ideal = y_p_pri / h_p_ref
    g_cal = calibration_retention * g_ideal
    y_d_guarded = y_d_pri - g_cal * h_d_ref
    y_p_guarded = y_p_pri - g_cal * h_p_ref
    p_d_guarded = float(abs(y_d_guarded) ** 2)
    p_p_guarded = float(abs(y_p_guarded) ** 2)
    lcr_guarded = 10.0 * math.log10(p_d_guarded / p_p_guarded)
    norm2_guarded = norm2_hr_dsb + float(abs(g_cal) ** 2)
    wng_guarded = 10.0 * math.log10(p_d_guarded / norm2_guarded)
    dt_loss_guarded = -10.0 * math.log10(p_d_guarded / p_d_hr_dsb)
    dt_worst_phase_loss = -20.0 * math.log10(1.0 - abs(g_cal * h_d_ref) / abs(y_d_pri))
    verdict_guarded = upper_limit_verdict(
        [[-lcr_guarded, -lcr_guarded], [-wng_guarded, -wng_guarded], [dt_loss_guarded, dt_loss_guarded]],
        upper_limits,
    )

    # 2d: Unguarded cross-zone canceller (unfrozen scalar Wiener adaptation during speech)
    g_unfrozen_equal_dt = (
        y_d_pri * h_d_ref.conjugate() + y_p_pri * h_p_ref.conjugate()
    ) / (abs(h_d_ref) ** 2 + abs(h_p_ref) ** 2 + sensor_noise_power)
    y_d_unfrozen_dt = y_d_pri - g_unfrozen_equal_dt * h_d_ref
    dt_loss_unguarded_equal_db = -10.0 * math.log10(float(abs(y_d_unfrozen_dt) ** 2) / p_d_hr_dsb)

    g_unfrozen_driver_only = (y_d_pri * h_d_ref.conjugate()) / (abs(h_d_ref) ** 2 + sensor_noise_power)
    y_d_unfrozen_driver = y_d_pri - g_unfrozen_driver_only * h_d_ref
    driver_only_self_cancel_db = -10.0 * math.log10(float(abs(y_d_unfrozen_driver) ** 2) / p_d_hr_dsb)

    norm2_unguarded = norm2_hr_dsb + float(abs(g_unfrozen_equal_dt) ** 2)
    wng_unguarded = 10.0 * math.log10(float(abs(y_d_unfrozen_dt) ** 2) / norm2_unguarded)
    verdict_unguarded = upper_limit_verdict(
        [
            [-lcr_guarded, -lcr_guarded],
            [-wng_unguarded, -wng_unguarded],
            [dt_loss_unguarded_equal_db, dt_loss_unguarded_equal_db],
        ],
        upper_limits,
    )

    candidates: dict[str, dict[str, object]] = {
        "ceiling_dsb": {
            "label": "顶棚中央双麦 DSB",
            "positions_m": mics_ceiling.tolist(),
            "driver_distances_m": r_d_ceil.tolist(),
            "passenger_distances_m": r_p_ceil.tolist(),
            "weights_real_imag": _complex_pairs(w_ceil),
            "weight_norm_squared": norm2_ceil,
            "driver_output_power": p_d_ceil,
            "passenger_leakage_power": p_p_ceil,
            "LCR_dB": lcr_ceil,
            "WNG_dB": wng_ceil,
            "dt_target_loss_db": 0.0,
            "verdict": verdict_ceil,
        },
        "headrest_dsb": {
            "label": "主驾头枕双麦 DSB",
            "positions_m": mics_hr1.tolist(),
            "driver_distances_m": r_d_hr1.tolist(),
            "passenger_distances_m": r_p_hr1.tolist(),
            "weights_real_imag": _complex_pairs(w_hr_dsb),
            "weight_norm_squared": norm2_hr_dsb,
            "driver_output_power": p_d_hr_dsb,
            "passenger_leakage_power": p_p_hr_dsb,
            "LCR_dB": lcr_hr_dsb,
            "WNG_dB": wng_hr_dsb,
            "dt_target_loss_db": 0.0,
            "verdict": verdict_hr_dsb,
        },
        "headrest_lcmv": {
            "label": "主驾头枕双麦零陷 LCMV",
            "positions_m": mics_hr1.tolist(),
            "weights_real_imag": _complex_pairs(w_hr_lcmv),
            "weight_norm_squared": norm2_hr_lcmv,
            "driver_output_power": p_d_hr_lcmv,
            "ideal_point_source_passenger_power": p_p_hr_lcmv_ideal,
            "passenger_leakage_power": p_p_hr_lcmv_pert,
            "LCR_dB": lcr_hr_lcmv,
            "WNG_dB": wng_hr_lcmv,
            "dt_target_loss_db": 0.0,
            "verdict": verdict_hr_lcmv,
        },
        "headrest_dsb_cross_cancel_unguarded": {
            "label": "头枕 DSB + 未加双讲保护跨区抵消",
            "weight_norm_squared": norm2_unguarded,
            "LCR_dB": lcr_guarded,
            "WNG_dB": wng_unguarded,
            "dt_target_loss_db": dt_loss_unguarded_equal_db,
            "driver_only_self_cancel_loss_db": driver_only_self_cancel_db,
            "verdict": verdict_unguarded,
        },
        "headrest_dsb_cross_cancel_guarded": {
            "label": "头枕 DSB + A2B 跨区参考抵消（含功率比冻结保护）",
            "reference_mic_position_m": mics_hr2[0].tolist(),
            "reference_driver_distance_m": float(r_d_hr2[0]),
            "reference_passenger_distance_m": float(r_p_hr2[0]),
            "calibrated_cancel_gain_abs": float(abs(g_cal)),
            "calibrated_cancel_gain_real_imag": _complex_pairs(g_cal),
            "weight_norm_squared": norm2_guarded,
            "driver_output_power": p_d_guarded,
            "passenger_leakage_power": p_p_guarded,
            "LCR_dB": lcr_guarded,
            "WNG_dB": wng_guarded,
            "dt_target_loss_db": dt_loss_guarded,
            "dt_worst_phase_loss_db": dt_worst_phase_loss,
            "verdict": verdict_guarded,
        },
    }

    eligible = [cid for cid, row in candidates.items() if row["verdict"] == "pass"]
    if len(eligible) != 1:
        raise RuntimeError(f"expected unique eligible cabin architecture, got {eligible}")

    return {
        "parameters": {
            "frequency_hz": frequency_hz,
            "sound_speed_m_s": sound_speed,
            "nominal_mouth_distance_m": nominal_mouth_m,
            "driver_position_m": s_driver.tolist(),
            "passenger_position_m": s_passenger.tolist(),
            "sensor_noise_power": sensor_noise_power,
            "cabin_scattering_perturbation": scattering_pert,
            "cross_zone_calibration_retention": calibration_retention,
            "hard_limits": {
                "LCR_min_dB": lcr_min_db,
                "WNG_min_dB": wng_min_db,
                "dt_target_loss_max_dB": dt_loss_max_db,
            },
        },
        "bus_clocks": bus_clocks,
        "bus_verdicts": bus_verdicts,
        "candidates": candidates,
        "eligible_architectures": eligible,
        "selected_architecture": eligible[0],
        "scope": "single-frequency spherical point-source + 5% perturbation control; not full cabin RIR",
    }


def icc_comb_and_howling_metrics(
    dsp_delay_ms: float,
    *,
    spk_to_ear_delay_ms: float = 1.0,
    direct_delay_ms: float = 5.0,
    reinforcement_ratio_g: float = 0.8,
    band_hz: tuple[float, float] = (300.0, 3400.0),
    feedback_coupling_db: float = -18.0,
    min_gain_margin_db: float = 6.0,
    p1150_limit_ms: float = 15.0,
) -> dict[str, object]:
    """Compute exact ICC excess delay, comb notches in band, and howling margin."""
    t_dsp_f = finite_real_scalar(dsp_delay_ms, "dsp_delay_ms")
    t_spk_f = finite_real_scalar(spk_to_ear_delay_ms, "spk_to_ear_delay_ms")
    t_dir_f = finite_real_scalar(direct_delay_ms, "direct_delay_ms")
    g = finite_real_scalar(reinforcement_ratio_g, "reinforcement_ratio_g")
    beta_db = finite_real_scalar(feedback_coupling_db, "feedback_coupling_db")
    gm_min_db = finite_real_scalar(min_gain_margin_db, "min_gain_margin_db")
    limit_ms = finite_real_scalar(p1150_limit_ms, "p1150_limit_ms")

    if t_dsp_f < 0.0 or t_spk_f < 0.0 or t_dir_f <= 0.0:
        raise ValueError("delays must be non-negative and direct_delay_ms must be positive")
    if not (0.0 < g < 1.0):
        raise ValueError("reinforcement_ratio_g must lie strictly in (0, 1)")
    if beta_db >= 0.0 or gm_min_db < 0.0 or limit_ms <= 0.0:
        raise ValueError("feedback_coupling_db must be < 0, min_gain_margin_db >= 0, p1150_limit_ms > 0")
    if not isinstance(band_hz, Sequence) or isinstance(band_hz, (str, bytes)) or len(band_hz) != 2:
        raise ValueError("band_hz must be (f_min_hz, f_max_hz)")
    f_min = _exact_decimal_fraction(band_hz[0], "band_hz[0]")
    f_max = _exact_decimal_fraction(band_hz[1], "band_hz[1]")
    if f_min <= 0 or f_min >= f_max:
        raise ValueError("band_hz must satisfy 0 < f_min < f_max")

    delta_rx_ms_frac = (
        _exact_decimal_fraction(t_dsp_f, "dsp_delay_ms")
        + _exact_decimal_fraction(t_spk_f, "spk_to_ear_delay_ms")
        - _exact_decimal_fraction(t_dir_f, "direct_delay_ms")
    )
    if delta_rx_ms_frac <= 0:
        raise ValueError("excess reinforcement delay Delta_Rx must be strictly positive")

    delta_rx_s_frac = delta_rx_ms_frac / Fraction(1000)
    spacing_hz_frac = Fraction(1) / delta_rx_s_frac
    half_spacing_frac = spacing_hz_frac / Fraction(2)

    # Destructive comb notches occur at f_m = (2*m + 1) / (2 * Delta_Rx), m >= 0.
    min_odd = math.ceil(f_min / half_spacing_frac)
    if min_odd % 2 == 0:
        min_odd += 1
    min_odd = max(1, min_odd)
    max_odd = math.floor(f_max / half_spacing_frac)
    if max_odd % 2 == 0:
        max_odd -= 1

    notches_frac = [half_spacing_frac * odd for odd in range(min_odd, max_odd + 1, 2)]
    p1150_pass = delta_rx_ms_frac < _exact_decimal_fraction(limit_ms, "p1150_limit_ms")
    ripple_db = 20.0 * math.log10((1.0 + g) / (1.0 - g))
    max_safe_gain_db = -beta_db - gm_min_db

    return {
        "dsp_delay_ms": float(t_dsp_f),
        "total_reinforced_arrival_ms": float(t_dsp_f + t_spk_f),
        "excess_delay_ms": float(delta_rx_ms_frac),
        "excess_delay_s_fraction": str(delta_rx_s_frac),
        "comb_spacing_hz": float(spacing_hz_frac),
        "comb_spacing_hz_fraction": str(spacing_hz_frac),
        "notch_count": len(notches_frac),
        "first_notch_hz": float(notches_frac[0]) if notches_frac else None,
        "last_notch_hz": float(notches_frac[-1]) if notches_frac else None,
        "notch_frequencies_hz": [float(v) for v in notches_frac],
        "p1150_pass": p1150_pass,
        "p1150_verdict": "pass" if p1150_pass else "fail",
        "peak_power_gain": (1.0 + g) ** 2,
        "notch_power_gain": (1.0 - g) ** 2,
        "peak_to_notch_ripple_db": ripple_db,
        "howling_max_gain_db": max_safe_gain_db,
        "howling_max_gain_linear": 10.0 ** (max_safe_gain_db / 20.0),
    }


def icc_and_certification_selection_case() -> dict[str, object]:
    """Return the E11-29 ITU-T P.1150 ICC and Teams/3QUEST certification control."""
    path_a = icc_comb_and_howling_metrics(12.0)
    path_b = icc_comb_and_howling_metrics(26.0)

    feedback_coupling_db = -18.0
    min_gain_margin_db = 6.0
    gain_checks: dict[str, dict[str, float | str]] = {}
    for gain_db in (12.0, 15.0, 18.0):
        loop_gain_db = gain_db + feedback_coupling_db
        margin_db = -loop_gain_db
        if loop_gain_db >= 0.0:
            regime = "unstable_howling"
            verdict = "fail"
        elif margin_db < min_gain_margin_db:
            regime = "insufficient_margin_ringing"
            verdict = "fail"
        else:
            regime = "safe_gain_margin"
            verdict = "pass"
        gain_checks[f"{gain_db:.1f}dB"] = {
            "icc_gain_db": gain_db,
            "loop_gain_db": loop_gain_db,
            "gain_margin_db": margin_db,
            "regime": regime,
            "verdict": verdict,
        }

    # Enterprise conferencing speakerphone certification (Teams / 3QUEST / TCLw):
    # Convert lower-bound requirements x >= L into minimized intervals [-x_high, -x_low] <= -L.
    cert_metrics = {
        "3quest_s_mos": {"measured": [4.1, 4.1], "min_limit": 3.8},
        "3quest_n_mos": {"measured": [4.3, 4.3], "min_limit": 3.8},
        "3quest_g_mos": {"measured": [4.0, 4.0], "min_limit": 3.8},
        "tclw_nominal_db": {"measured": [58.0, 58.0], "min_limit": 55.0},
        "tclw_max_volume_db": {"measured": [44.5, 47.0], "min_limit": 46.0},
        "send_dsp_latency_ms": {"measured": [28.0, 32.0], "max_limit": 40.0},
    }
    minimized_intervals: list[list[float]] = []
    minimized_limits: list[float] = []
    per_metric_verdicts: dict[str, str] = {}
    for key, spec in cert_metrics.items():
        lo, hi = spec["measured"]
        if "min_limit" in spec:
            iv = [-hi, -lo]
            lim = -float(spec["min_limit"])
        else:
            iv = [lo, hi]
            lim = float(spec["max_limit"])
        minimized_intervals.append(iv)
        minimized_limits.append(lim)
        per_metric_verdicts[key] = upper_limit_verdict([iv], [lim])

    teams_verdict = upper_limit_verdict(minimized_intervals, minimized_limits)

    return {
        "parameters": {
            "direct_delay_ms": 5.0,
            "spk_to_ear_delay_ms": 1.0,
            "reinforcement_ratio_g": 0.8,
            "speech_band_hz": [300.0, 3400.0],
            "p1150_excess_delay_limit_ms": 15.0,
            "feedback_coupling_db": feedback_coupling_db,
            "min_gain_margin_db": min_gain_margin_db,
        },
        "path_A_8ms": path_a,
        "path_B_22ms": path_b,
        "notch_density_ratio_B_over_A": path_b["notch_count"] / path_a["notch_count"],
        "peak_to_notch_ripple_db": path_a["peak_to_notch_ripple_db"],
        "howling_max_gain_db": path_a["howling_max_gain_db"],
        "howling_max_gain_linear": path_a["howling_max_gain_linear"],
        "howling_gain_candidates": gain_checks,
        "teams_certification": {
            "metrics": cert_metrics,
            "per_metric_verdicts": per_metric_verdicts,
        },
        "teams_certification_verdict": teams_verdict,
        "scope": "closed-form two-ray comb ripple, scalar Barkhausen margin, and interval certification table",
    }
