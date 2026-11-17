"""Deterministic 16 kHz synthetic wearable-array audio fixtures and 3-tier scoring.

Provides six synchronized 32000-sample signals sharing ``common_export_gain = 1.0``
and steady scoring window ``[1600, 30400)`` (28800 samples = 1.8 s, exact integer
cycles of all target and disturbance tones). Analytic steady-state expectations,
double-precision float64 measurements, and quantized PCM16 integer measurements
are strictly separated.
"""
from __future__ import annotations

import numpy as np

from codes.chapters.ch00.core.audio_samples import delay_samples
from codes.chapters.ch02.core.conventions import validate_waveforms
from codes.chapters.ch16.core.wearable_fusion import cross_modal_mvdr_mmse

SAMPLE_RATE = 16000
SAMPLES = 32000
FADE_SAMPLES = 1600
SCORE_START = 1600
SCORE_STOP = 30400
SCORING_SAMPLES = SCORE_STOP - SCORE_START
HEARTHROUGH_DELAY_SAMPLES = 32  # Exact 2.0 ms at 16 kHz -> 250 Hz notch, 4000 Hz peak
COMMON_EXPORT_GAIN = 1.0

TARGET_LOW_HZ = 250
TARGET_HIGH_HZ = 4000
TARGET_AMP = 0.10  # Each tone steady power = 0.005; combined clean power = 0.010

FILE_NAMES = {
    "wearable_clean_ref": "wearable_clean_ref.wav",
    "wearable_air_wind_2ch": "wearable_air_wind_2ch.wav",
    "wearable_inear_occluded": "wearable_inear_occluded.wav",
    "wearable_inear_equalized": "wearable_inear_equalized.wav",
    "wearable_fused_output": "wearable_fused_output.wav",
    "wearable_hearthrough_comb": "wearable_hearthrough_comb.wav",
}

CHANNEL_COUNTS = {
    "wearable_clean_ref": 1,
    "wearable_air_wind_2ch": 2,
    "wearable_inear_occluded": 1,
    "wearable_inear_equalized": 1,
    "wearable_fused_output": 1,
    "wearable_hearthrough_comb": 1,
}

LIMITS = (
    "Deterministic dual-tone target (250 Hz and 4000 Hz) with finite-window orthogonal "
    "disturbance tones, known occlusion gains (2.0x at 250 Hz, 0.1x at 4000 Hz), and "
    "exact 32-sample (2 ms) causal delay; not speech, real wind turbulence, measured "
    "ear-canal transfer functions, blind RTF estimation, or a subjective listening study. "
    "No per-file normalization, clipping, dither, or post-hoc gain/delay fitting."
)


def complex_pairs(value) -> list:
    arr = np.asarray(value, dtype=np.complex128)
    if not np.all(np.isfinite(arr)):
        raise ValueError("complex metadata must be finite")
    return np.stack([arr.real, arr.imag], axis=-1).tolist()


def raised_cosine_envelope(n_samples: int = SAMPLES, fade_samples: int = FADE_SAMPLES) -> np.ndarray:
    """Smooth C^1 raised-cosine edge fade, identically 1.0 on [fade_samples, n_samples - fade_samples)."""
    if (
        isinstance(n_samples, (bool, np.bool_))
        or not isinstance(n_samples, (int, np.integer))
        or isinstance(fade_samples, (bool, np.bool_))
        or not isinstance(fade_samples, (int, np.integer))
        or n_samples < 1
        or fade_samples < 0
        or 2 * fade_samples > n_samples
    ):
        raise ValueError("n_samples >= 1 and 0 <= 2 * fade_samples <= n_samples required")
    env = np.ones(int(n_samples), dtype=np.float64)
    if fade_samples == 0:
        return env
    head = np.arange(fade_samples, dtype=np.float64)
    env[:fade_samples] = 0.5 * (1.0 - np.cos(np.pi * head / fade_samples))
    tail_dist = (n_samples - 1) - np.arange(n_samples - fade_samples, n_samples, dtype=np.float64)
    env[n_samples - fade_samples :] = 0.5 * (1.0 - np.cos(np.pi * tail_dist / fade_samples))
    return env


def parameters() -> dict:
    low_sol = cross_modal_mvdr_mmse(2.0, np.diag([4.0, 0.04]), phi_s=1.0)
    high_sol = cross_modal_mvdr_mmse(0.1, np.diag([0.05, 0.04]), phi_s=1.0)
    return {
        "exercise_ids": ["E16-12", "E16-13", "E16-14", "E16-16"],
        "sample_rate_hz": SAMPLE_RATE,
        "samples_per_channel": SAMPLES,
        "duration_s": 2.0,
        "fade_samples_each_end": FADE_SAMPLES,
        "envelope": "0.5*(1-cos(pi*r/1600)) on [0,1600) and [30400,32000), 1.0 on [1600,30400)",
        "scoring_interval_samples": [SCORE_START, SCORE_STOP],
        "scoring_samples_per_channel": SCORING_SAMPLES,
        "target_frequencies_hz": [TARGET_LOW_HZ, TARGET_HIGH_HZ],
        "target_amplitudes": [TARGET_AMP, TARGET_AMP],
        "target_band_steady_powers": [0.005, 0.005],
        "target_total_steady_power": 0.01,
        "air_wind_low_frequencies_hz": [125, 375],
        "air_wind_low_amplitudes": [0.20, 0.20],
        "air_high_noise_frequencies_hz": [3500, 4500],
        "air_high_noise_quadrature_amplitudes": [0.02, 0.01],
        "inear_occlusion_gains": {"low_250hz": 2.0, "high_4000hz": 0.1},
        "inear_sensor_noise_frequencies_hz": {"low": 175, "high": 5000},
        "inear_sensor_noise_amplitudes": {"low": 0.02, "high": 0.02},
        "naive_inverse_equalization_gains": {"low": 0.5, "high": 10.0},
        "mvdr_subband_weights": {
            "low_air_inear": [float(np.real(v)) for v in low_sol["w_mvdr"]],
            "high_air_inear": [float(np.real(v)) for v in high_sol["w_mvdr"]],
        },
        "hearthrough_delay_samples": HEARTHROUGH_DELAY_SAMPLES,
        "hearthrough_delay_ms": 2.0,
        "hearthrough_gains": [1.0, 1.0],
        "common_export_gain": COMMON_EXPORT_GAIN,
    }


def build_signals_and_components() -> tuple[dict[str, np.ndarray], dict[str, dict[str, np.ndarray]]]:
    """Construct exact 6 CxN float64 waveforms and subband target/noise components."""
    n = np.arange(SAMPLES, dtype=np.float64)
    env = raised_cosine_envelope(SAMPLES, FADE_SAMPLES)

    s_low = env * (TARGET_AMP * np.cos(2.0 * np.pi * TARGET_LOW_HZ * n / SAMPLE_RATE))
    s_high = env * (TARGET_AMP * np.cos(2.0 * np.pi * TARGET_HIGH_HZ * n / SAMPLE_RATE))
    s_clean = s_low + s_high

    # Outer 2-mic air observations: wind-like low-freq orthogonal tones + slight high-freq tones.
    w_air0_low = env * (0.20 * np.cos(2.0 * np.pi * 125.0 * n / SAMPLE_RATE))
    w_air1_low = env * (0.20 * np.cos(2.0 * np.pi * 375.0 * n / SAMPLE_RATE))
    v_air0_high = env * (
        0.02 * np.cos(2.0 * np.pi * 3500.0 * n / SAMPLE_RATE)
        + 0.01 * np.sin(2.0 * np.pi * 3500.0 * n / SAMPLE_RATE)
    )
    v_air1_high = env * (
        0.02 * np.cos(2.0 * np.pi * 4500.0 * n / SAMPLE_RATE)
        + 0.01 * np.sin(2.0 * np.pi * 4500.0 * n / SAMPLE_RATE)
    )

    air_ch0 = s_clean + w_air0_low + v_air0_high
    air_ch1 = s_clean + w_air1_low + v_air1_high
    air_2ch = np.vstack([air_ch0, air_ch1])

    # Occluded in-ear / bone observation: 2.0x boost at 250 Hz, 0.1x attenuation at 4000 Hz.
    v_in_low = env * (0.02 * np.cos(2.0 * np.pi * 175.0 * n / SAMPLE_RATE))
    v_in_high = env * (0.02 * np.cos(2.0 * np.pi * 5000.0 * n / SAMPLE_RATE))
    inear_low = 2.0 * s_low + v_in_low
    inear_high = 0.1 * s_high + v_in_high
    inear_occ = inear_low + inear_high

    # Naive inverse equalization: 0.5x on low band, 10.0x on high band (+20 dB high-freq noise!).
    inear_eq_clean = 0.5 * (2.0 * s_low) + 10.0 * (0.1 * s_high)
    inear_eq_noise = 0.5 * v_in_low + 10.0 * v_in_high
    inear_eq = 0.5 * inear_low + 10.0 * inear_high

    # Frequency-selective cross-modal MVDR fusion combining air_ch0 and inear_occ:
    # Low band: [1/401, 200/401]^T; High band: [80/81, 10/81]^T.
    g_low = cross_modal_mvdr_mmse(2.0, np.diag([4.0, 0.04]), phi_s=1.0)["w_mvdr"].real
    g_high = cross_modal_mvdr_mmse(0.1, np.diag([0.05, 0.04]), phi_s=1.0)["w_mvdr"].real
    air0_low = s_low + w_air0_low
    air0_high = s_high + v_air0_high

    fused_low = g_low[0] * air0_low + g_low[1] * inear_low
    fused_high = g_high[0] * air0_high + g_high[1] * inear_high
    fused = fused_low + fused_high
    fused_clean = (g_low[0] + 2.0 * g_low[1]) * s_low + (g_high[0] + 0.1 * g_high[1]) * s_high
    fused_noise = (g_low[0] * w_air0_low + g_low[1] * v_in_low) + (
        g_high[0] * v_air0_high + g_high[1] * v_in_high
    )

    # Hear-through comb filtering: direct leak + 32-sample (2.0 ms) delayed path.
    hearthrough = s_clean + delay_samples(s_clean, HEARTHROUGH_DELAY_SAMPLES)

    signals = {
        "wearable_clean_ref": s_clean[np.newaxis, :],
        "wearable_air_wind_2ch": air_2ch,
        "wearable_inear_occluded": inear_occ[np.newaxis, :],
        "wearable_inear_equalized": inear_eq[np.newaxis, :],
        "wearable_fused_output": fused[np.newaxis, :],
        "wearable_hearthrough_comb": hearthrough[np.newaxis, :],
    }
    components = {
        "wearable_clean_ref": {"clean": s_clean, "noise": np.zeros_like(s_clean)},
        "wearable_air_wind_2ch": {"clean": s_clean, "noise": w_air0_low + v_air0_high},
        "wearable_inear_occluded": {"clean": 2.0 * s_low + 0.1 * s_high, "noise": v_in_low + v_in_high},
        "wearable_inear_equalized": {"clean": inear_eq_clean, "noise": inear_eq_noise},
        "wearable_fused_output": {"clean": fused_clean, "noise": fused_noise},
        "wearable_hearthrough_comb": {
            "clean": hearthrough,
            "noise": np.zeros_like(s_clean),
        },
    }
    return signals, components


def analytic_expectations() -> dict:
    """Closed-form steady-state powers, tone amplitudes, and reference MSEs."""
    low_noise = 0.005 / 100.25  # 0.020/(401^2) + 40000*0.0002/(401^2)
    high_noise = 0.005 / 20.25  # 6400*0.00025/(81^2) + 100*0.0002/(81^2)
    fused_mse = low_noise + high_noise
    # Occluded clean target vs clean_ref: (2-1)^2 * 0.005 + (0.1-1)^2 * 0.005 = 0.005 + 0.00405 = 0.00905
    occ_dist = 0.005 + 0.81 * 0.005
    occ_noise = 0.0002 + 0.0002
    # Equalized noise: 0.25 * 0.0002 + 100 * 0.0002 = 0.00005 + 0.020 = 0.02005
    eq_noise = 0.00005 + 0.020
    # Hear-through steady state: 250 Hz -> 0 (mse 0.005), 4000 Hz -> 2*s_high (power 0.020, dist (2-1)^2*0.005 = 0.005)
    return {
        "wearable_clean_ref": {
            "steady_mean_square_ch0": 0.010,
            "tone_complex_amplitudes_ch0": [[0.10, 0.0], [0.10, 0.0]],
            "target_distortion_power": 0.0,
            "residual_noise_power": 0.0,
            "reference_mse_ch0": 0.0,
        },
        "wearable_air_wind_2ch": {
            "steady_mean_square_ch0": 0.010 + 0.020 + 0.00025,
            "tone_complex_amplitudes_ch0": [[0.10, 0.0], [0.10, 0.0]],
            "target_distortion_power": 0.0,
            "residual_noise_power": 0.02025,
            "reference_mse_ch0": 0.02025,
        },
        "wearable_inear_occluded": {
            "steady_mean_square_ch0": 0.020 + 0.00005 + occ_noise,
            "tone_complex_amplitudes_ch0": [[0.20, 0.0], [0.01, 0.0]],
            "target_distortion_power": occ_dist,
            "residual_noise_power": occ_noise,
            "reference_mse_ch0": occ_dist + occ_noise,
        },
        "wearable_inear_equalized": {
            "steady_mean_square_ch0": 0.010 + eq_noise,
            "tone_complex_amplitudes_ch0": [[0.10, 0.0], [0.10, 0.0]],
            "target_distortion_power": 0.0,
            "residual_noise_power": eq_noise,
            "reference_mse_ch0": eq_noise,
        },
        "wearable_fused_output": {
            "steady_mean_square_ch0": 0.010 + fused_mse,
            "tone_complex_amplitudes_ch0": [[0.10, 0.0], [0.10, 0.0]],
            "target_distortion_power": 0.0,
            "residual_noise_power": fused_mse,
            "reference_mse_ch0": fused_mse,
        },
        "wearable_hearthrough_comb": {
            "steady_mean_square_ch0": 0.020,
            "tone_complex_amplitudes_ch0": [[0.0, 0.0], [0.20, 0.0]],
            "target_distortion_power": 0.010,
            "residual_noise_power": 0.0,
            "reference_mse_ch0": 0.010,
        },
    }


def _tone_projections(channel_slice: np.ndarray) -> np.ndarray:
    idx = np.arange(SCORE_START, SCORE_STOP, dtype=np.float64)
    freqs = np.array([TARGET_LOW_HZ, TARGET_HIGH_HZ], dtype=np.float64)
    phasors = np.exp(-2.0j * np.pi * freqs[:, np.newaxis] * idx[np.newaxis, :] / SAMPLE_RATE)
    return (2.0 / SCORING_SAMPLES) * (phasors @ channel_slice)


def measure_signal(
    samples: np.ndarray,
    key: str,
    *,
    reference: np.ndarray,
    components: dict[str, np.ndarray] | None = None,
    pcm: bool = False,
) -> dict:
    """Score float64 or decoded PCM16 waveform on [1600, 30400)."""
    if type(pcm) is not bool:
        raise ValueError("pcm must be bool")
    if key not in FILE_NAMES:
        raise ValueError(f"unknown wearable audio key: {key}")
    x = validate_waveforms(samples)
    ref = validate_waveforms(reference)
    expected_ch = CHANNEL_COUNTS[key]
    if x.shape != (expected_ch, SAMPLES) or ref.shape != (1, SAMPLES):
        raise ValueError(f"unexpected waveform shape for {key}")

    win = x[:, SCORE_START:SCORE_STOP]
    ref_win = ref[0, SCORE_START:SCORE_STOP]
    ms_by_ch = [float(np.mean(ch**2)) for ch in win]
    proj_ch0 = _tone_projections(win[0])
    diff_ch0 = win[0] - ref_win
    mse_ch0 = float(np.mean(diff_ch0**2))
    ref_power = float(np.mean(ref_win**2))

    result = {
        "scoring_samples_per_channel": SCORING_SAMPLES,
        "peak_all_samples": float(np.max(np.abs(x))),
        "mean_square_by_channel": ms_by_ch,
        "tone_projections_ch0_real_imag": complex_pairs(proj_ch0),
        "reference_mse_ch0": mse_ch0,
        "reference_nmse_ch0": float(mse_ch0 / ref_power),
    }

    if not pcm and components is not None:
        clean_win = components["clean"][SCORE_START:SCORE_STOP]
        noise_win = components["noise"][SCORE_START:SCORE_STOP]
        dist_win = clean_win - ref_win
        result["float_decomposition_ch0"] = {
            "target_distortion_power": float(np.mean(dist_win**2)),
            "residual_noise_power": float(np.mean(noise_win**2)),
            "cross_term": float(2.0 * np.mean(dist_win * noise_win)),
            "reconstruction_max_abs_error": float(np.max(np.abs(win[0] - clean_win - noise_win))),
        }

    if pcm:
        x_int = np.rint(x * 32768.0).astype(np.int64)
        ref_int = np.rint(ref[0] * 32768.0).astype(np.int64)
        win_int = x_int[:, SCORE_START:SCORE_STOP]
        ref_win_int = ref_int[SCORE_START:SCORE_STOP]
        sq_sums = [int(np.sum(ch * ch)) for ch in win_int]
        err_int = win_int[0] - ref_win_int
        err_sq_sum = int(np.sum(err_int * err_int))
        ref_sq_sum = int(np.sum(ref_win_int * ref_win_int))
        result["pcm_integer_metrics"] = {
            "pcm_decode_divisor": 32768,
            "integer_power_sample_denominator": SCORING_SAMPLES,
            "integer_power_scale_denominator": 32768**2,
            "integer_scored_squared_sum_by_channel": sq_sums,
            "integer_reference_error_squared_sum_ch0": err_sq_sum,
            "integer_reference_squared_sum": ref_sq_sum,
            "integer_nmse_ch0": float(err_sq_sum / ref_sq_sum),
        }

    return result


def run_experiment() -> tuple[dict, dict[str, np.ndarray]]:
    signals, components = build_signals_and_components()
    ref = signals["wearable_clean_ref"]
    float_scores = {
        key: measure_signal(arr, key, reference=ref, components=components[key], pcm=False)
        for key, arr in signals.items()
    }
    report = {
        "parameters": parameters(),
        "analytic_expectations": analytic_expectations(),
        "float_measurements": float_scores,
        "limits": LIMITS,
    }
    return report, signals
