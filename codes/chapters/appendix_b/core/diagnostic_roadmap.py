"""Cascade bottleneck attribution and troubleshooting roadmap kernels for Appendix B.

Implements E17-16:
1. Exact SINR / dB accounting across a 3-disturbance cascade (Linear AEC ->
   MIMO-WPE -> Spatial MVDR), contrasting downstream-only tuning under an
   unfixed upstream echo bottleneck against upstream-first remediation.
2. Deterministic upstream-to-downstream fault triage mapping 6 canonical
   telemetry signatures across Chapters 6, 7, 8, 10, 13, 14, and 15.

Pure Python/NumPy implementation with strict finite non-bool input validation.
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


def _stage_record(
    name: str,
    target_power: float,
    echo_power: float,
    reverb_power: float,
    sensor_noise_power: float,
    baseline_sinr_db: float | None = None,
    prev_sinr_db: float | None = None,
) -> dict[str, object]:
    disturbance = echo_power + reverb_power + sensor_noise_power
    if disturbance <= 0.0 or not math.isfinite(disturbance):
        raise ValueError(f"{name} total disturbance power must be positive and finite")
    sinr_linear = target_power / disturbance
    sinr_db = 10.0 * math.log10(sinr_linear)
    components = {
        "echo": echo_power,
        "late_reverb": reverb_power,
        "sensor_noise": sensor_noise_power,
    }
    dominant = max(components.items(), key=lambda item: item[1])[0]
    record: dict[str, object] = {
        "stage": name,
        "target_power": target_power,
        "echo_power": echo_power,
        "late_reverb_power": reverb_power,
        "sensor_noise_power": sensor_noise_power,
        "total_disturbance_power": disturbance,
        "sinr_linear": sinr_linear,
        "sinr_db": sinr_db,
        "dominant_bottleneck": dominant,
        "dominant_share": components[dominant] / disturbance,
    }
    if prev_sinr_db is not None:
        record["step_gain_db"] = sinr_db - prev_sinr_db
    if baseline_sinr_db is not None:
        record["cumulative_gain_db"] = sinr_db - baseline_sinr_db
    return record


def cascade_bottleneck_budget(
    *,
    target_power: float = 1.0,
    echo_power: float = 10.0,
    late_reverb_power: float = 0.50,
    sensor_noise_power: float = 0.01,
    aec_erle_db: float = 20.0,
    wpe_reverb_retention: float = 0.20,
    wrong_mvdr_diffuse_factor: float = 0.50,
    full_mvdr_factor: float = 0.25,
) -> dict[str, object]:
    """Compute exact 3-stage cascade SINR attribution (E17-16)."""
    p_s = _validate_real_scalar(target_power, "target_power")
    p_echo = _validate_real_scalar(echo_power, "echo_power")
    p_rev = _validate_real_scalar(late_reverb_power, "late_reverb_power")
    p_v = _validate_real_scalar(sensor_noise_power, "sensor_noise_power")
    erle_db = _validate_real_scalar(aec_erle_db, "aec_erle_db")
    wpe_ret = _validate_real_scalar(wpe_reverb_retention, "wpe_reverb_retention")
    wrong_bf = _validate_real_scalar(
        wrong_mvdr_diffuse_factor, "wrong_mvdr_diffuse_factor"
    )
    full_bf = _validate_real_scalar(full_mvdr_factor, "full_mvdr_factor")

    if p_s <= 0.0 or p_echo < 0.0 or p_rev < 0.0 or p_v <= 0.0:
        raise ValueError("target/sensor powers must be > 0 and echo/reverb >= 0")
    if erle_db < 0.0:
        raise ValueError("aec_erle_db must be >= 0")
    if not (0.0 < wpe_ret <= 1.0):
        raise ValueError("wpe_reverb_retention must lie in (0, 1]")
    if not (0.0 < wrong_bf <= 1.0) or not (0.0 < full_bf <= 1.0):
        raise ValueError("MVDR power factors must lie in (0, 1]")

    raw = _stage_record("raw_single_mic", p_s, p_echo, p_rev, p_v)
    raw_db = float(raw["sinr_db"])

    # Wrong triage: AEC broken (ERLE = 0 dB) and loudspeaker near target direction,
    # while downstream WPE + 2-mic MVDR suppress late reverb and sensor noise.
    wrong_reverb = p_rev * wpe_ret * wrong_bf
    wrong_noise = p_v * wrong_bf
    wrong = _stage_record(
        "wrong_downstream_only_triage",
        p_s,
        p_echo,
        wrong_reverb,
        wrong_noise,
        baseline_sinr_db=raw_db,
        prev_sinr_db=raw_db,
    )

    # Correct upstream-first cascade:
    aec_scale = 10.0 ** (-erle_db / 10.0)
    echo_after_aec = p_echo * aec_scale
    stage1_aec = _stage_record(
        "after_linear_aec",
        p_s,
        echo_after_aec,
        p_rev,
        p_v,
        baseline_sinr_db=raw_db,
        prev_sinr_db=raw_db,
    )
    reverb_after_wpe = p_rev * wpe_ret
    stage2_wpe = _stage_record(
        "after_mimo_wpe",
        p_s,
        echo_after_aec,
        reverb_after_wpe,
        p_v,
        baseline_sinr_db=raw_db,
        prev_sinr_db=float(stage1_aec["sinr_db"]),
    )
    stage3_mvdr = _stage_record(
        "after_4mic_mvdr",
        p_s,
        echo_after_aec * full_bf,
        reverb_after_wpe * full_bf,
        p_v * full_bf,
        baseline_sinr_db=raw_db,
        prev_sinr_db=float(stage2_wpe["sinr_db"]),
    )

    return {
        "inputs": {
            "target_power": p_s,
            "echo_power": p_echo,
            "late_reverb_power": p_rev,
            "sensor_noise_power": p_v,
            "raw_ser_db": 10.0 * math.log10(p_s / p_echo) if p_echo > 0.0 else math.inf,
            "raw_drr_db": 10.0 * math.log10(p_s / p_rev) if p_rev > 0.0 else math.inf,
            "raw_snr_db": 10.0 * math.log10(p_s / p_v),
            "aec_erle_db": erle_db,
            "wpe_reverb_retention": wpe_ret,
            "wrong_mvdr_diffuse_factor": wrong_bf,
            "full_mvdr_factor": full_bf,
        },
        "raw": raw,
        "wrong_downstream_only": wrong,
        "upstream_first_cascade": [stage1_aec, stage2_wpe, stage3_mvdr],
        "wrong_vs_correct_final_gap_db": float(stage3_mvdr["sinr_db"])
        - float(wrong["sinr_db"]),
    }


_DEFAULT_TELEMETRY: dict[str, float] = {
    "sro_ppm": 0.0,
    "hear_through_latency_us": 50.0,
    "low_freq_wind_coherence": 0.90,
    "low_freq_level_db": -35.0,
    "far_end_active_ratio": 0.0,
    "far_end_erle_db": 25.0,
    "aec_ref_delay_error_ms": 0.0,
    "rt60_s": 0.25,
    "drr_db": 6.0,
    "wpe_prediction_delay_frames": 3.0,
    "overlap_speech_ratio": 0.02,
    "der_percent": 6.0,
}


def diagnose_pipeline_symptom(telemetry: dict[str, object]) -> dict[str, object]:
    """Deterministically triage the highest-priority upstream fault from telemetry."""
    if not isinstance(telemetry, dict):
        raise ValueError("telemetry must be a dict")
    unknown = set(telemetry) - set(_DEFAULT_TELEMETRY)
    if unknown:
        raise ValueError(f"unsupported telemetry keys: {sorted(unknown)}")

    merged = dict(_DEFAULT_TELEMETRY)
    for key, raw_val in telemetry.items():
        merged[key] = _validate_real_scalar(raw_val, f"telemetry[{key!r}]")

    if merged["hear_through_latency_us"] < 0.0:
        raise ValueError("telemetry['hear_through_latency_us'] must be non-negative")
    if merged["rt60_s"] < 0.0:
        raise ValueError("telemetry['rt60_s'] must be non-negative")
    if merged["wpe_prediction_delay_frames"] < 0.0:
        raise ValueError("telemetry['wpe_prediction_delay_frames'] must be non-negative")
    if merged["der_percent"] < 0.0:
        raise ValueError("telemetry['der_percent'] must be non-negative")
    for ratio_key in ("low_freq_wind_coherence", "far_end_active_ratio", "overlap_speech_ratio"):
        if not 0.0 <= merged[ratio_key] <= 1.0:
            raise ValueError(f"telemetry[{ratio_key!r}] must lie in [0, 1]")

    # Priority 1: Hardware sample clock / cross-device SRO synchronization
    if abs(merged["sro_ppm"]) > 5.0:
        return {
            "fault_id": "clock_sro_drift",
            "priority_rank": 1,
            "stage": "clock_and_bus_sync",
            "chapter_pointer": "§10.2 / §13.15",
            "action": (
                "Lock shared hardware clock (TDM/PDM/PTPv2) or run continuous SRO "
                "estimation and polyphase resampling before covariance estimation"
            ),
        }

    # Priority 2: Hear-through electro-acoustic comb-filtering latency budget
    if merged["hear_through_latency_us"] > 100.0:
        comb_first_notch_hz = 1.0 / (2.0 * merged["hear_through_latency_us"] * 1e-6)
        return {
            "fault_id": "hear_through_latency_comb",
            "priority_rank": 2,
            "stage": "hear_through_codec_path",
            "chapter_pointer": "§14.10",
            "first_comb_notch_hz": comb_first_notch_hz,
            "action": (
                "Reduce converter + DSP hear-through path to <= 100 us or attenuate "
                "active playback in the passive-leakage overlap band"
            ),
        }

    # Priority 3: Outdoor wearable wind buffeting / incoherent low-frequency overload
    if merged["low_freq_wind_coherence"] < 0.25 and merged["low_freq_level_db"] > -15.0:
        return {
            "fault_id": "wearable_wind_saturation",
            "priority_rank": 3,
            "stage": "sensor_mechanical_and_cross_modal",
            "chapter_pointer": "§14.8 / §14.9",
            "action": (
                "Detect wind via low spatial coherence, enforce WNG / high-pass "
                "protection, and fuse occluded in-ear microphone below crossover"
            ),
        }

    # Priority 4: Linear AEC reference alignment / ERLE before dereverb or beamforming
    if (
        merged["far_end_active_ratio"] > 0.10 and merged["far_end_erle_db"] < 10.0
    ) or abs(merged["aec_ref_delay_error_ms"]) > 2.0:
        return {
            "fault_id": "aec_reference_or_erle",
            "priority_rank": 4,
            "stage": "linear_aec",
            "chapter_pointer": "§6.1 / §6.8",
            "action": (
                "Align loopback reference bulk delay, verify system post-processing "
                "tap point, and check speaker clipping before tuning WPE or MVDR"
            ),
        }

    # Priority 5: MIMO-WPE late reverberation & guard-delay configuration
    if merged["rt60_s"] > 0.40 and (
        merged["drr_db"] < 0.0 or merged["wpe_prediction_delay_frames"] < 2.0
    ):
        return {
            "fault_id": "late_reverberation_wpe",
            "priority_rank": 5,
            "stage": "mimo_wpe_dereverberation",
            "chapter_pointer": "§7.2 / §7.6",
            "action": (
                "Configure WPE guard delay Delta >= 2 frames to protect direct/early "
                "speech reflections and match filter tap count K to room RT60"
            ),
        }

    # Priority 6: Multi-speaker overlap diarization & GSS target leakage
    if merged["overlap_speech_ratio"] > 0.10 and merged["der_percent"] > 20.0:
        return {
            "fault_id": "overlap_diarization_gss",
            "priority_rank": 6,
            "stage": "separation_and_diarization",
            "chapter_pointer": "§8.4 / §15.5",
            "action": (
                "Upgrade single-label clustering to multi-label spatial TS-VAD and "
                "drive GSS covariance estimation with exclusive target masks"
            ),
        }

    return {
        "fault_id": "nominal",
        "priority_rank": 0,
        "stage": "healthy",
        "chapter_pointer": "§11.2 / §17.4",
        "action": "All upstream physical preconditions satisfied; proceed to worst-case scenario acceptance",
    }


def diagnostic_roadmap_controls() -> dict[str, object]:
    """Return JSON-safe E17-16 cascade attribution and fault-routing fixture."""
    budget = cascade_bottleneck_budget()
    canonical_cases = {
        "clock_drift": diagnose_pipeline_symptom({"sro_ppm": 100.0}),
        "hear_through_comb": diagnose_pipeline_symptom(
            {"hear_through_latency_us": 2000.0}
        ),
        "wind_buffeting": diagnose_pipeline_symptom(
            {"low_freq_wind_coherence": 0.08, "low_freq_level_db": -6.0}
        ),
        "broken_aec_reference": diagnose_pipeline_symptom(
            {"far_end_active_ratio": 0.6, "far_end_erle_db": 1.5, "aec_ref_delay_error_ms": 25.0}
        ),
        "wpe_overwhitening_or_tail": diagnose_pipeline_symptom(
            {"rt60_s": 0.75, "drr_db": -3.0, "wpe_prediction_delay_frames": 1.0}
        ),
        "meeting_overlap_der": diagnose_pipeline_symptom(
            {"overlap_speech_ratio": 0.24, "der_percent": 28.5}
        ),
        "compound_clock_plus_overlap": diagnose_pipeline_symptom(
            {"sro_ppm": 80.0, "overlap_speech_ratio": 0.24, "der_percent": 32.0}
        ),
        "healthy_pipeline": diagnose_pipeline_symptom({}),
    }
    return {
        "exercise_id": "E17-16",
        "cascade_budget": budget,
        "symptom_triage": canonical_cases,
        "scope": (
            "deterministic power-budget attribution and upstream-first diagnostic "
            "triage across Ch 6/7/8/10/13/14/15; not a blind estimator or live telemetry daemon"
        ),
    }
