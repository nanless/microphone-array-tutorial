"""Run AEC, WPE, separation, direction tracking and beamforming on one recording."""

from __future__ import annotations

import argparse
import hashlib
import json

import numpy as np

from codes.chapters.ch02.core.spectral import istft, stft
from codes.chapters.ch00.core.audio_samples import pcm16_bytes, read_pcm16
from codes.chapters.ch03.core.geometry import plane_wave_steering
from codes.chapters.ch04.core.doa import srp_phat
from codes.chapters.ch05.core.beamforming import dsb_weights
from codes.chapters.ch06.core.aec import nlms
from codes.chapters.ch07.core.dereverberation import offline_wpe
from codes.chapters.ch08.core.gss_teaching import guided_cacgmm_mvdr
from codes.chapters.ch08.core.separation import pit_permutation, si_sdr
from codes.chapters.ch09.core.tracking import ConstantVelocityKalman, white_acceleration_covariance


RATE = 8000
SEED = 603
SAMPLES = 6 * RATE
FFT_LENGTH = 256
HOP = 64
MICROPHONES = np.array([[0., 0.], [.06, 0.]])
TRACK_FRAMES = 16


def _delay(signal, delay_samples):
    samples = np.arange(signal.size)
    return np.interp(samples - delay_samples, samples, signal, left=0., right=0.)


def _recording():
    rng = np.random.default_rng(SEED)
    carriers = rng.standard_normal((3, SAMPLES))
    # Distinct onset and amplitude patterns give the two sources different envelopes.
    envelopes = np.repeat([[0., 0., .5, .1, .8, .6, 1., .35, .6, .8, .4, .9],
                           [0., 0., 0., .9, .8, 1., 0., .9, .25, 1., .8, 0.]],
                          RATE // 2, axis=1)
    sources = .075 * envelopes * carriers[:2]
    # A fractional advance of less than one sample must not put the onset
    # inside the first microphone-only echo training interval.
    sources[:, RATE] = 0.
    target_angle = np.interp(np.arange(SAMPLES) / RATE, [1., 5.], [18., 42.])
    angles = np.stack((target_angle, np.full(SAMPLES, -32.)))
    sample_delay = -.06 * RATE / 343 * np.sin(np.deg2rad(angles))
    direct = np.stack((sources, np.array([_delay(s, d) for s, d in zip(sources, sample_delay)])), axis=1)
    reverberant = direct.sum(axis=0)
    for delay, gain in ((256, .60), (448, -.35), (640, .25)):
        reflection = np.stack((_delay(sources[0], delay) + .8 * _delay(sources[1], delay + 12),
                               .85 * _delay(sources[0], delay + 1.3) + _delay(sources[1], delay + 11.2)))
        reverberant += gain * reflection
    playback = .12 * carriers[2]
    echo = np.stack((.8 * _delay(playback, 4) - .28 * _delay(playback, 17),
                     .7 * _delay(playback, 5) + .22 * _delay(playback, 20)))
    observed = reverberant + echo + .002 * rng.standard_normal((2, SAMPLES))
    wav = pcm16_bytes(observed, RATE)
    return wav, playback, echo, direct, envelopes, target_angle


def _synthesis(spectrum):
    return istft(spectrum, n_fft=FFT_LENGTH, hop_length=HOP, length=SAMPLES)


def _track(spectrum, frequencies):
    sector = np.arange(0., 61.)  # Follow the talker on the positive side of the array.
    observations = []
    centers = []
    edges = range(RATE // HOP + TRACK_FRAMES, spectrum.shape[-1] - TRACK_FRAMES, TRACK_FRAMES)
    for start in edges:
        frame = spectrum[:, :, start:start + TRACK_FRAMES]
        score = srp_phat(frame, frequencies, MICROPHONES, np.deg2rad(sector))
        observations.append(sector[score.argmax()])
        centers.append(start * HOP + (TRACK_FRAMES - 1) * HOP // 2)

    interval = TRACK_FRAMES * HOP / RATE
    kalman = ConstantVelocityKalman([observations[0], 0.], np.diag([25., 100.]),
                                    white_acceleration_covariance(interval, 80.))
    estimates = [kalman.state[0]]
    for angle in observations[1:]:
        kalman.predict(interval)
        estimates.append(kalman.update(angle, 9.)[0])
    return np.array(centers, dtype=int), np.array(observations), np.array(estimates)


def _beamform(spectrum, frequencies, centers, angles):
    per_frame = np.interp(np.arange(spectrum.shape[-1]) * HOP, centers, angles)
    steer = plane_wave_steering(MICROPHONES, frequencies, np.deg2rad(per_frame))
    weights = dsb_weights(steer.reshape(-1, len(MICROPHONES))).reshape(steer.shape)
    output = np.einsum("cfj,jfc->fj", spectrum, weights.conj())
    return _synthesis(output[None, :, :])[0]


def run_experiment(*, blind=False):
    wav, playback, echo, direct, envelopes, true_angles = _recording()
    _, observed = read_pcm16(wav)
    freeze = np.arange(SAMPLES) >= RATE
    cancelled = [nlms(playback, channel, 32, freeze=freeze) for channel in observed]
    residual = np.stack([result[0] for result in cancelled])
    echo_hat = np.stack([result[1] for result in cancelled])
    score = slice(2 * RATE, 11 * RATE // 2)
    direct_at_mic_one = direct[:, 0, score]
    removed_echo = 10 * np.log10(np.sum(echo[:, score]**2) / np.sum((echo - echo_hat)[:, score]**2))
    spectrum = stft(residual, n_fft=FFT_LENGTH, hop_length=HOP)
    wpe = offline_wpe(spectrum.transpose(1, 0, 2), delay=3, taps=8, iterations=2).transpose(1, 0, 2)
    frequencies = np.fft.rfftfreq(FFT_LENGTH, 1 / RATE)
    wpe_waveform = _synthesis(wpe)

    samples_at_frame_center = np.minimum(np.arange(spectrum.shape[-1]) * HOP, SAMPLES - 1)
    activity = (envelopes[:, samples_at_frame_center] > 0).T
    guided = []
    for speaker in range(2):
        estimate = guided_cacgmm_mvdr(wpe.transpose(1, 0, 2), activity, target_speaker=speaker)
        guided.append(estimate["output"])
    guided_waveform = _synthesis(np.stack(guided))

    doa_results = {}
    for name, spatial in (("AEC", spectrum), ("AEC_WPE", wpe)):
        centers, measured, tracks = _track(spatial, frequencies)
        reference = true_angles[centers]
        beamformed = _beamform(spatial, frequencies, centers, tracks)
        doa_results[name] = {
            "observations": len(centers),
            "first_and_last_observation_center_samples": [int(centers[0]), int(centers[-1])],
            "angle_rmse_degrees_before_tracking": float(np.sqrt(np.mean((measured - reference)**2))),
            "angle_rmse_degrees_after_tracking": float(np.sqrt(np.mean((tracks - reference)**2))),
            "tracked_dsb_target_si_sdr_db": si_sdr(beamformed[score], direct_at_mic_one[0]),
        }

    def mono_scores(signal):
        return [si_sdr(signal[score], reference) for reference in direct_at_mic_one]

    result = {
        "seed": SEED, "sample_rate_hz": RATE, "microphone_samples": SAMPLES,
        "recording_wav_sha256": hashlib.sha256(wav).hexdigest(),
        "first_second_for_echo_training": True,
        "waveform_sample_interval": [score.start, score.stop],
        "known_echo_component_attenuation_db": float(removed_echo),
        "si_sdr_db_per_reference": {
            "microphone_one": mono_scores(observed[0]),
            "AEC_microphone_one": mono_scores(residual[0]),
            "AEC_WPE_microphone_one": mono_scores(wpe_waveform[0]),
            "AEC_WPE_guided_separation": [si_sdr(signal[score], reference) for signal, reference
                                         in zip(guided_waveform, direct_at_mic_one)],
        },
        "doa_and_beamforming": doa_results,
    }
    if blind:
        from codes.chapters.ch04.core.upstream_contracts import check_unchanged
        from codes.chapters.ch08.examples.audit_separation_upstream_interfaces import package_unchanged
        from codes.chapters.ch08.examples.reproduce_auxiva_reference import pinned_auxiva_class

        auxiva = pinned_auxiva_class()
        result["blind_AuxIVA"] = {}
        for name, observation in (("AEC", spectrum), ("AEC_WPE", wpe)):
            model = auxiva(spatial_algorithm="IP", scale_restoration=True, reference_id=0,
                           record_loss=False)
            separated = _synthesis(model(observation, n_iter=20))
            permutation, mean = pit_permutation(separated[:, score], direct_at_mic_one)
            result["blind_AuxIVA"][name] = {"output_to_reference": list(permutation),
                                          "mean_si_sdr_db": mean}
        result["blind_AuxIVA"]["source_commit"] = check_unchanged(pinned_auxiva_class.identity)["head"]
        package_unchanged(pinned_auxiva_class.installed)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blind", action="store_true", help="also use the pinned ssspy checkout")
    print(json.dumps(run_experiment(blind=parser.parse_args().blind), indent=2, allow_nan=False))
