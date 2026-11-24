"""Extend E11-20 to causal FIRs and measure their digital response."""

from __future__ import annotations

import json

import numpy as np

from codes.chapters.ch03.core.geometry import plane_wave_steering
from codes.chapters.ch05.core.beamforming import dsb_weights, mvdr_weights
from codes.chapters.ch11.core.selection_physics import beamformer_selection_case


def run_experiment() -> dict:
    case = beamformer_selection_case()["parameters"]
    positions = np.array(case["positions_m"])
    rate = 16_000
    fft_length = 4096
    delay = 64
    frequencies = np.fft.rfftfreq(fft_length, 1 / rate)
    band = (frequencies >= 300) & (frequencies <= 3400)

    def steering(key):
        return plane_wave_steering(
            positions, frequencies, np.deg2rad(case[key]),
            reference=case["reference_channel"], sound_speed=case["sound_speed_m_s"],
        )

    nominal = steering("nominal_target_deg")
    actual = steering("actual_target_deg")
    interferer = steering("interference_deg")
    noise = case["independent_sensor_noise_power"]
    covariance = interferer[:, :, None] * interferer[:, None, :].conj() + noise * np.eye(3)
    common_delay = np.exp(-2j * np.pi * frequencies * delay / rate)
    rows = []
    for loading in (None, 0., .1, 1.):
        weights = dsb_weights(nominal) if loading is None else mvdr_weights(
            covariance, nominal, relative_diagonal_loading=loading,
        )
        # The beamformer multiplies each input by the conjugate of its weight.
        elementary = np.fft.irfft(weights.conj(), n=fft_length, axis=0)
        taps = np.roll(elementary, delay, axis=0)[:2 * delay + 1]
        taps *= np.hamming(len(taps))[:, None]
        implemented = np.fft.rfft(taps, n=fft_length, axis=0)
        nominal_gain = np.sum(implemented * nominal, axis=1)
        actual_gain = np.sum(implemented * actual, axis=1) / common_delay
        interference_gain = np.sum(implemented * interferer, axis=1)
        norm = np.sum(np.abs(implemented)**2, axis=1)
        white_noise_gain = 10 * np.log10(np.abs(nominal_gain[band])**2 / norm[band])
        target_error = np.abs(actual_gain[band] - 1)
        output_noise = np.abs(interference_gain[band])**2 + noise * norm[band]
        measured_delay = int(np.argmax(np.abs(taps.sum(axis=1))))
        rows.append({
            "relative_loading": loading,
            "minimum_wng_db": float(white_noise_gain.min()),
            "maximum_target_error": float(target_error.max()),
            "maximum_target_error_frequency_hz": float(frequencies[band][np.argmax(target_error)]),
            "equal_bin_mean_nmse": float(np.mean(target_error**2 + output_noise)),
            "acoustic_verdict": (
                "pass" if white_noise_gain.min() >= case["WNG_lower_limit_dB"]
                and target_error.max() <= case["response_amplitude_error_upper_limit"] else "fail"
            ),
            "broadside_digital_impulse_peak_samples": measured_delay,
        })
    return {
        "sample_rate_hz": rate, "fft_length": fft_length,
        "taps_per_channel": 2 * delay + 1, "design_delay_samples": delay,
        "design_delay_ms": 1000 * delay / rate, "evaluated_band_hz": [300, 3400],
        "evaluated_bin_count": int(band.sum()), "candidates": rows,
        "device_end_to_end_latency_ms": None,
    }


if __name__ == "__main__":
    print(json.dumps(run_experiment(), indent=2, allow_nan=False))
