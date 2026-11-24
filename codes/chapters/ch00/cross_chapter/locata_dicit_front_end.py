"""Run five front-end stages on LOCATA development Task 2, DICIT recording 1."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from codes.chapters.ch00.core.audio_samples import pcm16_bytes
from codes.chapters.ch00.io_contracts import validate_asset_directory
from codes.chapters.ch02.core.spectral import istft, stft
from codes.chapters.ch03.core.geometry import plane_wave_steering
from codes.chapters.ch04.core.doa import srp_phat
from codes.chapters.ch04.core.upstream_contracts import check_unchanged
from codes.chapters.ch05.core.beamforming import dsb_weights
from codes.chapters.ch05.examples.locata_estimated_steering import MICROPHONES
from codes.chapters.ch06.core.aec import erle_db
from codes.chapters.ch07.core.dereverberation import offline_wpe
from codes.chapters.ch08.core.separation import pit_permutation, si_sdr
from codes.chapters.ch08.examples.audit_separation_upstream_interfaces import package_unchanged
from codes.chapters.ch08.examples.reproduce_auxiva_reference import pinned_auxiva_class
from codes.chapters.ch09.core.tracking import ConstantVelocityKalman, white_acceleration_covariance
from scripts.course_home import prepare_locata as locata


RATE = 16000
SAMPLES = 100843
N_FFT = 512
HOP = 128
SCAN_FRAMES = 8
GRID = np.arange(-90.0, 91.0)
FREQUENCIES = np.fft.rfftfreq(N_FFT, 1 / RATE)
BAND = (FREQUENCIES >= 375) & (FREQUENCIES <= 3000)
MIC = locata.REFERENCE
AUXIVA_PAIR = (MIC, 4)
LAGS = np.arange(-1600, 1601)
TAPS = 1024
LEADING_TAPS = 64
LOADING = 1e-3
A_FIT = (42400, 45600)  # 2.65–2.85 s
A_ONLY = (46560, 49600)  # 2.91–3.10 s
B_FIT = tuple(n // 3 for n in locata.CALIBRATION)
B_ONLY = tuple(n // 3 for n in locata.WINDOWS['B_only'])
OVERLAP = tuple(n // 3 for n in locata.WINDOWS['overlap'])
MEMBERS = {
    **locata.MEMBERS,
    'dicit/audio_source_loudspeaker3.wav': '607486982b2f87be3f7d3022432e14c05756164625c2b855ca05d6b143957526',
    'dicit/audio_source_loudspeaker4.wav': '9eb1addb5efb02b21d527e24b97e8f848f78a14af66465a00d1d4f03eb266cbf',
}


def _read_waveform(path: Path, channels: int) -> np.ndarray:
    data = subprocess.check_output([
        'ffmpeg', '-nostdin', '-v', 'error', '-i', str(path),
        '-ar', str(RATE), '-f', 'f64le', '-',
    ])
    waveform = np.frombuffer(data, dtype='<f8').reshape(-1, channels).T
    if waveform.shape != (channels, SAMPLES) or not np.isfinite(waveform).all():
        raise ValueError(f'Unexpected resampled waveform: {path}')
    return waveform


def _reference_delay(source: np.ndarray, mic: np.ndarray,
                     interval: tuple[int, int]) -> tuple[int, float]:
    """Fit one common integer shift on a single-source interval."""
    start, stop = interval
    copies = sliding_window_view(source, stop - start)[start - LAGS]
    target = mic[start:stop]
    denominator = np.linalg.norm(copies, axis=1) * np.linalg.norm(target)
    corr = np.divide(copies @ target, denominator, out=np.zeros_like(denominator), where=denominator > 0)
    peak = int(np.argmax(np.abs(corr)))
    return int(LAGS[peak]), float(corr[peak])


def _delay(signal: np.ndarray, samples: int) -> np.ndarray:
    output = np.zeros_like(signal)
    if samples >= 0:
        output[samples:] = signal[:len(signal) - samples]
    else:
        output[:samples] = signal[-samples:]
    return output


def _cancel(playback: np.ndarray, waveform: np.ndarray, lag: int) -> np.ndarray:
    """Fit one loaded causal FIR per microphone and replay with frozen weights."""
    # Leave room on either side of the strongest arrival. The array and
    # dry playback share a file timeline but the acoustic paths add delay.
    shifted = _delay(playback, lag - LEADING_TAPS)
    history = sliding_window_view(np.pad(shifted, (TAPS - 1, 0)), TAPS)
    fit = slice(*B_FIT)
    design = history[fit].copy()
    gram = design.T @ design
    gram[np.diag_indices(TAPS)] += LOADING * np.trace(gram) / TAPS
    # The row order in the design is oldest through newest.
    filters = np.linalg.solve(gram, design.T @ waveform[:, fit].T)[::-1].T
    fft_length = 1 << (SAMPLES + TAPS - 2).bit_length()
    echo = np.fft.irfft(np.fft.rfft(shifted, fft_length) * np.fft.rfft(filters, fft_length), fft_length)
    return waveform - echo[:, :SAMPLES]


def _source_angles(paths: dict[str, Path]) -> np.ndarray:
    array = np.median(np.loadtxt(paths['dicit/position_array_dicit.txt'], skiprows=1), axis=0)
    rotation = array[12:21].reshape(3, 3)
    sources = np.stack([
        np.median(np.loadtxt(paths[f'benchmark2/position_source_loudspeaker{i}.txt'], skiprows=1), axis=0)[6:9]
        for i in (3, 4)
    ])
    local = (sources - array[6:9]) @ rotation
    return np.degrees(np.arctan2(local[:, 0], local[:, 1]))


def _complete_frames(centers: np.ndarray, interval: tuple[int, int]) -> np.ndarray:
    return np.flatnonzero((centers - N_FFT // 2 >= interval[0]) & (centers + N_FFT // 2 <= interval[1]))


def _scan(spectrum: np.ndarray) -> np.ndarray:
    return srp_phat(spectrum, FREQUENCIES[BAND], MICROPHONES, np.deg2rad(GRID))


def _track(spectrum: np.ndarray, centers: np.ndarray, eligible: np.ndarray,
           initial_angle: float, target_mask: np.ndarray | None = None) -> list[dict]:
    support = spectrum[:, BAND]
    if target_mask is not None:
        support = np.where(target_mask[None, BAND] >= 0.5, support, 0)
    interval = SCAN_FRAMES * HOP / RATE
    process = white_acceleration_covariance(interval, 50.0)
    kalman = ConstantVelocityKalman([initial_angle, 0.0], np.diag([36.0, 64.0]), process)
    target_sector = np.abs(GRID - initial_angle) <= 22.0
    events = []
    for start in range(0, eligible.size - SCAN_FRAMES + 1, SCAN_FRAMES):
        frames = eligible[start:start + SCAN_FRAMES]
        global_peak = target_peak = None
        if np.any(support[:, :, frames]):
            score = _scan(support[:, :, frames])
            global_peak = float(GRID[np.argmax(score)])
            target_peak = float(GRID[target_sector][np.argmax(score[target_sector])])
        if events:
            kalman.predict(interval)
        if target_peak is not None:
            kalman.update(target_peak, 16.0)
        events.append({
            'state_sample': int((centers[frames[0]] + centers[frames[-1]]) // 2),
            'available_sample': int(centers[frames[-1]] + N_FFT // 2),
            'global_peak_deg': global_peak,
            'sector_peak_deg': target_peak,
            'filtered_deg': float(kalman.state[0]),
        })
    return events


def _synth(spectrum: np.ndarray) -> np.ndarray:
    return istft(spectrum, n_fft=N_FFT, hop_length=HOP, length=SAMPLES)


def _beam(spectrum: np.ndarray, centers: np.ndarray, initial_angle: float,
          track: Sequence[dict] = (), target_mask: np.ndarray | None = None) -> np.ndarray:
    edges = [event['available_sample'] for event in track]
    direction = np.array([initial_angle] + [event['filtered_deg'] for event in track])
    steering = plane_wave_steering(MICROPHONES, FREQUENCIES, np.deg2rad(direction), reference=MIC)
    weights = dsb_weights(steering.reshape(-1, len(MICROPHONES))).reshape(steering.shape)
    active = np.searchsorted(edges, centers, side='right')
    output = np.einsum('mft,tfm->ft', spectrum, weights[active].conj())
    if target_mask is not None:
        output *= target_mask
    return _synth(output[None])[0]


def _errors(events: list[dict], target: float) -> dict:
    summary = {}
    for key in ('global_peak_deg', 'sector_peak_deg', 'filtered_deg'):
        values = np.array([event[key] for event in events if event[key] is not None])
        summary[key] = {
            'observations': len(values),
            'rmse_deg': float(np.sqrt(np.mean((values - target) ** 2))) if values.size else None,
        }
    return summary


def run(input_dir: Path) -> tuple[dict, dict[str, np.ndarray]]:
    paths = {name: input_dir / locata.RECORDING / name for name in MEMBERS}
    for name, path in paths.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != MEMBERS[name]:
            raise ValueError(f'Unexpected LOCATA member: {name}')
    raw = _read_waveform(paths['dicit/audio_array_dicit.wav'], len(MICROPHONES))
    playback = [_read_waveform(paths[f'dicit/audio_source_loudspeaker{i}.wav'], 1)[0] for i in (3, 4)]
    a_active, b_active = [np.loadtxt(paths[f'dicit/VAD_dicit_loudspeaker{i}.txt'], skiprows=1).astype(bool)
                          for i in (3, 4)]
    for (start, stop), labels in ((A_FIT, a_active & ~b_active), (A_ONLY, a_active & ~b_active),
                                  (B_FIT, b_active & ~a_active), (B_ONLY, b_active & ~a_active),
                                  (OVERLAP, a_active & b_active)):
        if not labels[start * 3:stop * 3].all():
            raise ValueError('LOCATA VAD differs inside an experiment interval')

    delays = [_reference_delay(source, raw[MIC], interval) for source, interval in zip(playback, (A_FIT, B_FIT))]
    references = np.stack([_delay(source, pair[0]) for source, pair in zip(playback, delays)])
    residual = _cancel(playback[1], raw, delays[1][0])
    a_only, b_only, overlap = (slice(*iv) for iv in (A_ONLY, B_ONLY, OVERLAP))
    a_projected = references[0, a_only]
    retention = np.dot(residual[MIC, a_only], a_projected) / np.dot(raw[MIC, a_only], a_projected)

    raw_stft = stft(raw, n_fft=N_FFT, hop_length=HOP)
    aec_stft = stft(residual, n_fft=N_FFT, hop_length=HOP)

    def wpe(spectrum: np.ndarray) -> np.ndarray:
        result = offline_wpe(spectrum.transpose(1, 0, 2), delay=3, taps=5, iterations=2)
        return result.transpose(1, 0, 2)

    centers = np.arange(raw_stft.shape[-1]) * HOP
    a_fit_frames = _complete_frames(centers, A_FIT)
    overlap_frames = _complete_frames(centers, OVERLAP)
    seed = float(GRID[np.argmax(_scan(raw_stft[:, BAND][:, :, a_fit_frames]))])
    angles = _source_angles(paths)
    auxiva = pinned_auxiva_class()
    routes = {'raw': raw_stft, 'aec': aec_stft, 'wpe': wpe(raw_stft), 'aec_wpe': wpe(aec_stft)}

    def score(waveform: np.ndarray) -> list[float]:
        return [si_sdr(waveform[overlap], reference[overlap]) for reference in references]

    comparisons = {}
    wavs = {'raw_mic7.wav': raw[MIC], 'aec_mic7.wav': residual[MIC]}
    for name, spectrum in routes.items():
        model = auxiva(spatial_algorithm='IP', scale_restoration=True, reference_id=0, record_loss=False)
        separated = model(spectrum[list(AUXIVA_PAIR)], n_iter=20)
        separated_wavs = _synth(separated)
        permutation, mean_score = pit_permutation(separated_wavs[:, overlap], references[:, overlap])
        power = np.abs(separated) ** 2
        target_slot = int(np.argmax(power[:, :, a_fit_frames].sum(axis=(1, 2))))
        total = power.sum(axis=0)
        mask = np.divide(power[target_slot], total, out=np.zeros_like(total), where=total > 0)

        ordinary = _track(spectrum, centers, overlap_frames, seed)
        masked = _track(spectrum, centers, overlap_frames, seed, mask)
        mic = _synth(spectrum[[MIC]])[0]
        ds_fixed = _beam(spectrum, centers, seed)
        ds_tracked = _beam(spectrum, centers, seed, ordinary)
        ds_masked = _beam(spectrum, centers, seed, masked, mask)
        comparisons[name] = {
            'mic7_si_sdr_db_a_b': score(mic),
            'auxiva_output_to_reference': list(permutation),
            'auxiva_output_si_sdr_db_a_b': [score(waveform) for waveform in separated_wavs],
            'auxiva_pit_mean_si_sdr_db': mean_score,
            'auxiva_target_slot_from_a_fit': target_slot,
            'track': ordinary,
            'track_rmse': _errors(ordinary, angles[0]),
            'masked_track': masked,
            'masked_track_rmse': _errors(masked, angles[0]),
            'fixed_ds_si_sdr_db_a_b': score(ds_fixed),
            'tracked_ds_si_sdr_db_a_b': score(ds_tracked),
            'masked_tracked_ds_si_sdr_db_a_b': score(ds_masked),
        }
        if name == 'aec_wpe':
            wavs.update({'aec_wpe_mic7.wav': mic, 'aec_wpe_auxiva_a.wav': separated_wavs[target_slot],
                         'aec_wpe_tracked_ds.wav': ds_tracked, 'aec_wpe_auxiva_mask_tracked_ds.wav': ds_masked})

    source_commit = check_unchanged(pinned_auxiva_class.identity)['head']
    package_unchanged(pinned_auxiva_class.installed)
    report = {
        'dataset': 'https://doi.org/10.5281/zenodo.3630471',
        'dataset_license': 'ODC-By-1.0',
        'archive_md5_from_zenodo': 'd5a5417c3f6b2ed0e43581dd06f504e6',
        'recording': locata.RECORDING.rstrip('/'),
        'archive_members_sha256': {locata.RECORDING + name: sha for name, sha in MEMBERS.items()},
        'implementation_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'ffmpeg_version': subprocess.check_output(['ffmpeg', '-version'], text=True).splitlines()[0],
        'numpy_version': np.__version__,
        'sample_rate_hz': RATE, 'samples': SAMPLES, 'channels': raw.shape[0],
        'one_based_mic_reference': MIC + 1,
        'one_based_auxiva_pair': [i + 1 for i in AUXIVA_PAIR],
        'auxiva_source_commit': source_commit, 'auxiva_iterations': 20,
        'stft_fft': N_FFT, 'stft_hop': HOP, 'frames_per_observation': SCAN_FRAMES,
        'minimum_and_maximum_scan_hz': [375, 3000],
        'initial_a_deg': seed, 'optical_a_b_deg': angles.tolist(),
        'a_fit_samples': A_FIT, 'a_only_retention_samples': A_ONLY, 'b_fit_samples': B_FIT,
        'b_only_erle_samples': B_ONLY, 'double_talk_score_samples': OVERLAP,
        'double_talk_complete_stft_frames': len(overlap_frames),
        'playback_delays_samples_a_b': [pair[0] for pair in delays],
        'playback_alignment_correlation_a_b': [pair[1] for pair in delays],
        'diagnostic_held_out_a_b_delay_and_correlation': [
            _reference_delay(playback[0], raw[MIC], A_ONLY), _reference_delay(playback[1], raw[MIC], B_ONLY)],
        'aec_taps': TAPS, 'aec_leading_taps_relative_to_b_peak': LEADING_TAPS,
        'aec_relative_diagonal_loading': LOADING,
        'aec_b_only_input_sum_squares': np.sum(raw[:, b_only] ** 2, axis=1).tolist(),
        'aec_b_only_residual_sum_squares': np.sum(residual[:, b_only] ** 2, axis=1).tolist(),
        'aec_b_only_erle_db_per_mic': [erle_db(x[b_only], r[b_only]) for x, r in zip(raw, residual)],
        'aec_near_a_only_projected_amplitude_change_db': float(20 * np.log10(abs(retention))),
        'aec_double_talk_mix_sum_squares': float(np.sum(raw[:, overlap] ** 2)),
        'aec_double_talk_residual_sum_squares': float(np.sum(residual[:, overlap] ** 2)),
        'wpe': {'delay_frames': 3, 'taps': 5, 'iterations': 2},
        'comparison': comparisons,
    }
    return report, wavs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_dir', type=Path, help='directory containing the extracted LOCATA dev/ tree')
    parser.add_argument('--wav-dir', type=Path, help='optional new or empty folder for six mono examples')
    args = parser.parse_args()
    report, wavs = run(args.input_dir)
    if args.wav_dir is not None:
        directory = validate_asset_directory(args.wav_dir, wavs, check=False)
        if directory.exists() and any(directory.iterdir()):
            raise ValueError('Use a new or empty --wav-dir')
        # Prepare every clip before writing and use the same gain in all of them.
        gain = float(min(4.0, 0.98 / max(np.max(np.abs(samples)) for samples in wavs.values())))
        encoded = {name: pcm16_bytes(samples * gain, RATE) for name, samples in wavs.items()}
        directory.mkdir(parents=True, exist_ok=True)
        for name, contents in encoded.items():
            with (directory / name).open('xb') as output:
                output.write(contents)
        report['wav_gain'] = gain
        report['wav_sha256'] = {name: hashlib.sha256(contents).hexdigest() for name, contents in encoded.items()}
    print(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
