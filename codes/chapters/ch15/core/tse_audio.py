"""Deterministic mathematical near-field distance-bubble audio controls for Chapter 15.

Synthesizes six 16 kHz, 32000-sample (2.0 s) waveforms at common export gain 1.0:
1. ``tse_near_target_ref.wav`` (1 ch): 3200 Hz near-field target reference at (0.40 m, 0 deg).
2. ``tse_far_interferer_ref.wav`` (1 ch): 2500 Hz co-directional far-field interferer at (4.00 m, 0 deg).
3. ``tse_array_3ch.wav`` (3 ch): 3-mic ULA (d = 0.10 m) receiving spherical near target,
   co-directional far interferer, and orthogonal sensor floor tones (1100/1300/1700 Hz).
4. ``tse_farfield_dsb.wav`` (1 ch): conventional far-field broadside DSB [1/3, 1/3, 1/3] on d = 0.10 m.
5. ``tse_bubble_d10cm.wav`` (1 ch): near-field distance-bubble LCMV on d = 0.10 m.
6. ``tse_bubble_d2cm_illcond.wav`` (1 ch): regularized small-aperture (d = 0.02 m, reg = 2e-4)
   distance-bubble output illustrating WNG collapse and target attenuation when r_1 >> r_F.
"""

from __future__ import annotations

import numpy as np

from codes.chapters.ch02.core.conventions import SPEED_OF_SOUND, finite_real_array
from codes.chapters.ch15.core.distance_bubble import (
    broadside_distance_bubble_weights,
    broadside_two_tap_real_fir,
)

SAMPLE_RATE = 16000
SAMPLES = 32000
FADE_SAMPLES = 800
STEADY_START = 1600
STEADY_STOP = 30400
STEADY_SAMPLES = STEADY_STOP - STEADY_START  # 28800 samples = 1.8 s

TARGET_FREQ_HZ = 3200.0
INTERFERER_FREQ_HZ = 2500.0
SENSOR_FLOOR_FREQS_HZ = (1100.0, 1300.0, 1700.0)

TARGET_RANGE_M = 0.40
INTERFERER_RANGE_M = 4.00
SPACING_GOOD_M = 0.10
SPACING_TINY_M = 0.02
TINY_REG = 2e-4

TARGET_AMPLITUDE = 0.16
INTERFERER_AMPLITUDE = 0.16
SENSOR_FLOOR_AMPLITUDE = 0.012

FILE_NAMES = {
    "near_target_ref": "tse_near_target_ref.wav",
    "far_interferer_ref": "tse_far_interferer_ref.wav",
    "array_3ch": "tse_array_3ch.wav",
    "farfield_dsb": "tse_farfield_dsb.wav",
    "bubble_d10cm": "tse_bubble_d10cm.wav",
    "bubble_d2cm_illcond": "tse_bubble_d2cm_illcond.wav",
}

CHANNELS_BY_KEY = {
    "near_target_ref": 1,
    "far_interferer_ref": 1,
    "array_3ch": 3,
    "farfield_dsb": 1,
    "bubble_d10cm": 1,
    "bubble_d2cm_illcond": 1,
}

OUTPUT_KEYS = ("farfield_dsb", "bubble_d10cm", "bubble_d2cm_illcond")

LIMITS = (
    "Deterministic interior-frequency sinusoids (3200/2500/1100/1300/1700 Hz) with exact integer "
    "cycles on [1600, 30400) and known broadside spherical propagation; not natural speech, "
    "a measured room, blind source range estimation, or a subjective listening test. "
    "Digital amplitude squared is not Pa^2; common export gain 1.0 with no per-file normalization, "
    "fitted gain, or post-hoc delay alignment."
)


def parameters() -> dict:
    """Return exact JSON-serializable metadata describing the 6-file synthetic scene."""
    return {
        "exercise_ids": ["E15-03", "E15-05", "E15-12"],
        "sample_rate_hz": SAMPLE_RATE,
        "samples_per_channel": SAMPLES,
        "fade_samples_each_end": FADE_SAMPLES,
        "steady_scoring_window": [STEADY_START, STEADY_STOP],
        "steady_scoring_samples": STEADY_SAMPLES,
        "speed_of_sound_m_per_s": float(SPEED_OF_SOUND),
        "target_frequency_hz": TARGET_FREQ_HZ,
        "interferer_frequency_hz": INTERFERER_FREQ_HZ,
        "sensor_floor_frequencies_hz": [float(f) for f in SENSOR_FLOOR_FREQS_HZ],
        "target_range_m": TARGET_RANGE_M,
        "interferer_range_m": INTERFERER_RANGE_M,
        "broadside_angle_deg": 0.0,
        "ula_spacing_good_m": SPACING_GOOD_M,
        "ula_spacing_tiny_m": SPACING_TINY_M,
        "tiny_aperture_tikhonov_reg": TINY_REG,
        "target_peak_amplitude": TARGET_AMPLITUDE,
        "interferer_peak_amplitude": INTERFERER_AMPLITUDE,
        "sensor_floor_peak_amplitude": SENSOR_FLOOR_AMPLITUDE,
        "target_steady_mean_square": float(0.5 * TARGET_AMPLITUDE ** 2),
        "interferer_steady_mean_square": float(0.5 * INTERFERER_AMPLITUDE ** 2),
        "sensor_floor_steady_mean_square": float(0.5 * SENSOR_FLOOR_AMPLITUDE ** 2),
        "common_export_gain": 1.0,
    }


def _envelope() -> np.ndarray:
    idx = np.arange(SAMPLES, dtype=np.float64)
    return np.minimum(np.minimum(idx / FADE_SAMPLES, (SAMPLES - 1.0 - idx) / FADE_SAMPLES), 1.0)


def _apply_two_tap_filter(array_3ch: np.ndarray, taps: np.ndarray) -> np.ndarray:
    """Apply causal 2-tap real FIR ``h_m[0] x_m[n] + h_m[1] x_m[n - 1]`` summed across channels."""
    out = taps[:, 0] @ array_3ch
    out[1:] += taps[:, 1] @ array_3ch[:, :-1]
    return out[None, :]


def _analytic_model_expectations() -> dict:
    """Compute exact steady-window analytic component powers and MSE expectations."""
    p_tar = 0.5 * (TARGET_AMPLITUDE ** 2)
    p_int = 0.5 * (INTERFERER_AMPLITUDE ** 2)
    p_noi = 0.5 * (SENSOR_FLOOR_AMPLITUDE ** 2)

    def eval_taps(d_m: float, taps: np.ndarray) -> dict:
        xs = np.array([-d_m, 0.0, d_m], dtype=np.float64)
        k_tar = 2.0 * np.pi * TARGET_FREQ_HZ / SPEED_OF_SOUND
        k_int = 2.0 * np.pi * INTERFERER_FREQ_HZ / SPEED_OF_SOUND
        a_tar = np.exp(-1j * k_tar * (np.hypot(TARGET_RANGE_M, xs) - TARGET_RANGE_M))
        a_int = np.exp(-1j * k_int * (np.hypot(INTERFERER_RANGE_M, xs) - INTERFERER_RANGE_M))

        om_tar = 2.0 * np.pi * TARGET_FREQ_HZ / SAMPLE_RATE
        om_int = 2.0 * np.pi * INTERFERER_FREQ_HZ / SAMPLE_RATE
        h_tar = taps[:, 0] + taps[:, 1] * np.exp(-1j * om_tar)
        h_int = taps[:, 0] + taps[:, 1] * np.exp(-1j * om_int)

        resp_tar = complex(np.dot(h_tar, a_tar))
        resp_int = complex(np.dot(h_int, a_int))

        noi_gain_sq = 0.0
        for m, fn in enumerate(SENSOR_FLOOR_FREQS_HZ):
            om_n = 2.0 * np.pi * fn / SAMPLE_RATE
            hm = taps[m, 0] + taps[m, 1] * np.exp(-1j * om_n)
            noi_gain_sq += float(abs(hm) ** 2)

        out_tar_power = float(p_tar * (abs(resp_tar) ** 2))
        tar_dist_power = float(p_tar * (abs(resp_tar - 1.0) ** 2))
        out_int_power = float(p_int * (abs(resp_int) ** 2))
        out_noi_power = float(p_noi * noi_gain_sq)
        total_power = out_tar_power + out_int_power + out_noi_power
        total_mse = tar_dist_power + out_int_power + out_noi_power
        return {
            "mean_square_per_channel": [total_power],
            "target_output_power": out_tar_power,
            "target_distortion_power": tar_dist_power,
            "interferer_leakage_power": out_int_power,
            "sensor_noise_output_power": out_noi_power,
            "reference_mse": total_mse,
            "reference_nmse": float(total_mse / p_tar),
            "interferer_rejection_db": float(-10.0 * np.log10(max(abs(resp_int) ** 2, 1e-30))),
        }

    dsb_taps = np.column_stack((np.full(3, 1.0 / 3.0), np.zeros(3)))
    bub10 = broadside_distance_bubble_weights(SPACING_GOOD_M, TARGET_RANGE_M, TARGET_FREQ_HZ)
    taps10 = broadside_two_tap_real_fir(bub10["weights"], TARGET_FREQ_HZ, SAMPLE_RATE)
    bub02 = broadside_distance_bubble_weights(
        SPACING_TINY_M, TARGET_RANGE_M, TARGET_FREQ_HZ, reg=TINY_REG
    )
    taps02 = broadside_two_tap_real_fir(bub02["weights"], TARGET_FREQ_HZ, SAMPLE_RATE)

    arr_ch_power = float(p_tar + p_int + p_noi)
    return {
        "near_target_ref": {"mean_square_per_channel": [float(p_tar)]},
        "far_interferer_ref": {"mean_square_per_channel": [float(p_int)]},
        "array_3ch": {"mean_square_per_channel": [arr_ch_power, arr_ch_power, arr_ch_power]},
        "farfield_dsb": eval_taps(SPACING_GOOD_M, dsb_taps),
        "bubble_d10cm": eval_taps(SPACING_GOOD_M, taps10),
        "bubble_d2cm_illcond": eval_taps(SPACING_TINY_M, taps02),
    }


def measure_signal(
    waveform: np.ndarray,
    key: str,
    *,
    reference: np.ndarray | None = None,
    pcm: bool = False,
) -> dict:
    """Measure steady-window float64 or actual PCM16 statistics on ``[1600, 30400)``."""
    if type(pcm) is not bool:
        raise ValueError("pcm must be bool")
    if not isinstance(key, str) or key not in FILE_NAMES:
        raise ValueError("unknown tse_audio signal key")
    x = finite_real_array(waveform, "waveform")
    if x.ndim != 2 or x.shape[0] < 1 or x.shape[1] < STEADY_STOP:
        raise ValueError("waveform must have 2-D shape (channels, samples) covering steady window")
    ref_flat: np.ndarray | None = None
    if reference is not None:
        ref_arr = finite_real_array(reference, "reference")
        ref_flat = ref_arr.reshape(-1)
        if ref_flat.size < STEADY_STOP:
            raise ValueError("reference must cover steady window")
    win = x[:, STEADY_START:STEADY_STOP]
    t_win = np.arange(STEADY_START, STEADY_STOP, dtype=np.float64) / SAMPLE_RATE
    k_tar = np.exp(-1j * 2.0 * np.pi * TARGET_FREQ_HZ * t_win)
    k_int = np.exp(-1j * 2.0 * np.pi * INTERFERER_FREQ_HZ * t_win)

    phasor_tar = (2.0 / STEADY_SAMPLES) * (win @ k_tar)
    phasor_int = (2.0 / STEADY_SAMPLES) * (win @ k_int)
    mean_sq = [float(np.mean(row * row)) for row in win]

    result: dict = {
        "steady_window": [STEADY_START, STEADY_STOP],
        "scoring_samples_per_channel": STEADY_SAMPLES,
        "peak_all_samples": float(np.max(np.abs(x))),
        "mean_square_per_channel": mean_sq,
        "target_phasor_real_imag": [[float(z.real), float(z.imag)] for z in phasor_tar],
        "interferer_phasor_real_imag": [[float(z.real), float(z.imag)] for z in phasor_int],
    }

    if key in OUTPUT_KEYS and ref_flat is not None:
        ref_win = ref_flat[STEADY_START:STEADY_STOP]
        diff = win[0] - ref_win
        mse = float(np.mean(diff * diff))
        ref_pow = float(np.mean(ref_win * ref_win))
        if ref_pow <= 0.0:
            raise ValueError("positive reference energy required in scoring window")
        result["reference_mse"] = mse
        result["reference_nmse"] = float(mse / ref_pow)

    if pcm:
        ints = np.rint(x * 32768.0).astype(np.int64)[:, STEADY_START:STEADY_STOP]
        e_ch = [int(sum(int(v) * int(v) for v in row)) for row in ints]
        d_ch = int(STEADY_SAMPLES * (32768 ** 2))
        pcm_int: dict = {
            "pcm_decode_divisor": 32768,
            "integer_squared_sum_E_per_channel": e_ch,
            "integer_denominator_D_per_channel": d_ch,
            "integer_squared_sum_E_all_channels": int(sum(e_ch)),
            "integer_denominator_D_all_channels": int(x.shape[0] * d_ch),
            "mean_square_all_channels": float(sum(e_ch) / (x.shape[0] * d_ch)),
        }
        if key in OUTPUT_KEYS and ref_flat is not None:
            ref_ints = np.rint(ref_flat[STEADY_START:STEADY_STOP] * 32768.0).astype(np.int64)
            err_ints = ints[0] - ref_ints
            e_err = int(sum(int(v) * int(v) for v in err_ints))
            p_ref = int(sum(int(v) * int(v) for v in ref_ints))
            if p_ref <= 0:
                raise ValueError("positive quantized reference energy required in scoring window")
            pcm_int["integer_reference_error_E"] = e_err
            pcm_int["integer_reference_power_P"] = p_ref
            pcm_int["integer_reference_nmse"] = float(e_err / p_ref)
        result["pcm_integer_measurements"] = pcm_int

    return result


def run_experiment() -> tuple[dict, dict[str, np.ndarray]]:
    """Synthesize all six waveforms in float64 and return ``(report, arrays)``."""
    env = _envelope()
    t_sec = np.arange(SAMPLES, dtype=np.float64) / SAMPLE_RATE

    near_ref = (env * TARGET_AMPLITUDE * np.cos(2.0 * np.pi * TARGET_FREQ_HZ * t_sec))[None, :]
    far_ref = (env * INTERFERER_AMPLITUDE * np.cos(2.0 * np.pi * INTERFERER_FREQ_HZ * t_sec))[None, :]
    sensor_noises = np.vstack(
        [env * SENSOR_FLOOR_AMPLITUDE * np.sin(2.0 * np.pi * fn * t_sec) for fn in SENSOR_FLOOR_FREQS_HZ]
    )

    def build_array(d_m: float) -> np.ndarray:
        xs = np.array([-d_m, 0.0, d_m], dtype=np.float64)
        tau_tar = (np.hypot(TARGET_RANGE_M, xs) - TARGET_RANGE_M) / SPEED_OF_SOUND
        tau_int = (np.hypot(INTERFERER_RANGE_M, xs) - INTERFERER_RANGE_M) / SPEED_OF_SOUND
        sig_tar = env[None, :] * TARGET_AMPLITUDE * np.cos(
            2.0 * np.pi * TARGET_FREQ_HZ * (t_sec[None, :] - tau_tar[:, None])
        )
        sig_int = env[None, :] * INTERFERER_AMPLITUDE * np.cos(
            2.0 * np.pi * INTERFERER_FREQ_HZ * (t_sec[None, :] - tau_int[:, None])
        )
        return (sig_tar + sig_int + sensor_noises).astype(np.float64)

    arr_good = build_array(SPACING_GOOD_M)
    arr_tiny = build_array(SPACING_TINY_M)

    dsb = np.mean(arr_good, axis=0, keepdims=True)
    bub_good = broadside_distance_bubble_weights(SPACING_GOOD_M, TARGET_RANGE_M, TARGET_FREQ_HZ)
    taps_good = broadside_two_tap_real_fir(bub_good["weights"], TARGET_FREQ_HZ, SAMPLE_RATE)
    out_good = _apply_two_tap_filter(arr_good, taps_good)

    bub_tiny = broadside_distance_bubble_weights(
        SPACING_TINY_M, TARGET_RANGE_M, TARGET_FREQ_HZ, reg=TINY_REG
    )
    taps_tiny = broadside_two_tap_real_fir(bub_tiny["weights"], TARGET_FREQ_HZ, SAMPLE_RATE)
    out_tiny = _apply_two_tap_filter(arr_tiny, taps_tiny)

    arrays = {
        "near_target_ref": near_ref,
        "far_interferer_ref": far_ref,
        "array_3ch": arr_good,
        "farfield_dsb": dsb,
        "bubble_d10cm": out_good,
        "bubble_d2cm_illcond": out_tiny,
    }

    analytic = _analytic_model_expectations()
    float_scores = {
        key: measure_signal(arr, key, reference=near_ref[0], pcm=False)
        for key, arr in arrays.items()
    }
    report = {
        "parameters": parameters(),
        "analytic_expectations": analytic,
        "float_measurements": float_scores,
        "limits": LIMITS,
    }
    return report, arrays
