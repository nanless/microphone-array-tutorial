"""Estimate relative clock drift, resample, and recover after a marked packet gap.

The default run uses a synthetic signal with separately known clock parameters.
Pass the official LibriWASN ZIP or its extracted record directory to run the
same estimator on three independently recorded phone channels.
"""

from __future__ import annotations

import argparse
from collections import deque
import hashlib
import io
import json
from pathlib import Path
import wave
import zipfile

import numpy as np


FS = 16_000
BLOCK = 8_192
LOOKBACK = 12
OBSERVATIONS = 4
MAX_LAG = 128
MAX_PPM = 500
MAX_RESIDUAL = 2.0
BUFFER = MAX_LAG + 2
GAP = (12, 15)
REAL_FIRST_BLOCK = 94
REAL_BLOCKS = 25

_WAV_MEMBER_ROOT = (
    "LibriWASN/200/0L/overlap_ratio_0.0_sil2.9_3.0_session8_actual0.0/record/"
)
_DEVICES = ("Pixel7", "Pixel6a", "Xiaomi")
_WAV_SHA256 = {
    "Pixel7": "815a0a766f32831e698f9df18c8fb3756573a417ce55bd719d857caf658e1ed6",
    "Pixel6a": "bbe1d6cc9fb29dddc170b21aa6a1bebfd4369e3f0d2cf4cc773ef6ff70986dd7",
    "Xiaomi": "bd99e1f46608cb8b2c53060baabfb154639e28dce3374ca3a2c43e0f2a824336",
}
_WAV_BYTES = {"Pixel7": 20_114_030, "Pixel6a": 20_114_254, "Xiaomi": 20_105_716}


def _relative_lag(reference: np.ndarray, remote: np.ndarray) -> float | None:
    """Return the band-limited PHAT peak, positive when remote needs a later index."""
    n_fft = 2 * BLOCK
    window = np.hanning(BLOCK)
    frequencies = np.fft.rfftfreq(n_fft, 1 / FS)
    cross = (
        np.fft.rfft(remote * window, n_fft)
        * np.fft.rfft(reference * window, n_fft).conj()
    )
    magnitude = np.abs(cross)
    use = (frequencies >= 300) & (frequencies <= 4_000) & (magnitude > 0)
    phase = np.zeros_like(cross)
    phase[use] = cross[use] / magnitude[use]
    circular = np.fft.irfft(phase)
    lags = np.arange(-MAX_LAG, MAX_LAG + 1)
    candidates = circular[lags]
    peak = int(np.argmax(candidates))
    if peak == 0 or peak == candidates.size - 1:
        return None
    left, value, right = candidates[peak - 1:peak + 2]
    competing = np.abs(lags - lags[peak]) > 8
    if value < 0.12 or value < 1.1 * np.max(candidates[competing]):
        return None
    return float(lags[peak] + 0.5 * (left - right) / (left - 2 * value + right))


def _fit(observations: deque[tuple[int, float]]) -> tuple[float, float]:
    """Median pair slope and median intercept in sample-index coordinates."""
    times, lags = np.asarray(observations).T
    first, second = np.triu_indices(len(times), k=1)
    slope = float(np.median((lags[second] - lags[first]) / (times[second] - times[first])))
    intercept = float(np.median(lags - times * slope))
    return slope, intercept


def _correlation(x: np.ndarray, y: np.ndarray) -> float:
    a = x - np.mean(x)
    b = y - np.mean(y)
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))


def _score(target: np.ndarray, aligned: np.ndarray, fixed: np.ndarray,
           *, clean_reference: bool) -> dict:
    score = {
        "samples": int(target.size),
        "tracked_correlation": _correlation(target, aligned),
        "fixed_correlation": _correlation(target, fixed),
    }
    if clean_reference:
        power = np.dot(target, target)
        for label, estimate in (("tracked", aligned), ("fixed", fixed)):
            residual = estimate - target
            score[f"{label}_nmse"] = float(np.dot(residual, residual) / power)
    return score


def _read_positions(remote: np.ndarray, positions: np.ndarray,
                    latest: int, received: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    first = np.floor(positions).astype(int)
    second = first + 1
    inside = (first >= 0) & (second < latest)
    valid = np.zeros(positions.size, dtype=bool)
    valid[inside] = received[first[inside]] & received[second[inside]]
    fraction = positions[valid] - first[valid]
    return valid, (1 - fraction) * remote[first[valid]] + fraction * remote[second[valid]]


def _run_pair(reference: np.ndarray, remote: np.ndarray,
              truth: np.ndarray | None = None, *,
              received: np.ndarray | None = None,
              arrival_ns: np.ndarray | None = None) -> dict:
    observations = deque(maxlen=LOOKBACK)
    model = None
    fixed_lag = None
    blocks = []
    targets, corrected, fixed = [], [], []
    if received is None:
        received = np.ones(remote.size, dtype=bool)
        received[GAP[0] * BLOCK:GAP[1] * BLOCK] = False
    for number, start in enumerate(range(0, reference.size, BLOCK)):
        stop = start + BLOCK
        center = start + BLOCK // 2
        deadline = stop + BUFFER
        row = {"block": number, "start_sample": start, "available_sample": deadline}
        blocks.append(row)
        visible = received if arrival_ns is None else received & (arrival_ns <= deadline * 1_000_000_000 // FS)
        if not visible[start:stop].all():
            row["status"] = "missing"
            observations.clear()
            model = None
            continue
        delay = _relative_lag(reference[start:stop], remote[start:stop])
        row["observed_lag_samples"] = delay
        if delay is None:
            row["status"] = "no_peak"
            observations.clear()
            model = None
            continue
        if model is not None and abs(delay - (model[0] * center + model[1])) > MAX_RESIDUAL:
            row["status"] = "lag_jump"
            observations.clear()
            model = None
        elif model is None:
            row["status"] = "acquiring"
        else:
            indices = np.arange(start, stop)
            offset = model[0] * indices + model[1]
            latest = min(deadline, remote.size)
            valid, aligned = _read_positions(remote, indices + offset, latest, visible)
            fixed_valid, baseline = _read_positions(remote, indices + fixed_lag, latest, visible)
            common = valid & fixed_valid
            targets.append((reference if truth is None else truth)[indices[common]])
            corrected.append(aligned[fixed_valid[valid]])
            fixed.append(baseline[valid[fixed_valid]])
            row.update(status="aligned", used_ppm=model[0] * 1e6,
                       predicted_lag_samples=model[0] * center + model[1],
                       aligned_samples=int(np.count_nonzero(valid)),
                       common_score_samples=int(np.count_nonzero(common)))
        observations.append((center, delay))
        if len(observations) >= OBSERVATIONS:
            fitted = _fit(observations)
            if abs(fitted[0]) <= MAX_PPM * 1e-6:
                model = fitted
                row["next_ppm"] = fitted[0] * 1e6
                if fixed_lag is None:
                    fixed_lag = fitted[0] * center + fitted[1]
            else:
                model = None

    aligned_rows = [row for row in blocks if row["status"] == "aligned"]
    after_gap = next((row for row in aligned_rows if row["block"] >= GAP[1] and row["aligned_samples"] > 0), None)
    result = {
        "fixed_lag_samples": fixed_lag,
        "aligned_blocks": len(aligned_rows),
        "median_used_ppm": float(np.median([row["used_ppm"] for row in aligned_rows]))
        if aligned_rows else None,
        "first_aligned_block_after_gap": after_gap["block"] if after_gap else None,
        "after_gap_first_aligned_output_s": (after_gap["available_sample"] - GAP[1] * BLOCK) / FS
        if after_gap else None,
        "tracked_center_residual_rms_samples": float(np.sqrt(np.mean([
            (row["predicted_lag_samples"] - row["observed_lag_samples"]) ** 2
            for row in aligned_rows
        ]))) if aligned_rows else None,
        "fixed_center_residual_rms_samples": float(np.sqrt(np.mean([
            (fixed_lag - row["observed_lag_samples"]) ** 2 for row in aligned_rows
        ]))) if aligned_rows else None,
        "score": _score(np.concatenate(targets), np.concatenate(corrected), np.concatenate(fixed),
                        clean_reference=truth is not None)
        if targets else None,
        "blocks": blocks,
    }
    return result


def _synthetic_signals():
    rng = np.random.default_rng(13022)
    length = REAL_BLOCKS * BLOCK
    index = np.arange(length, dtype=float)
    oversample = 16
    grid = np.arange(-64 * oversample, (length + 64) * oversample) / oversample
    frequencies = np.fft.rfftfreq(grid.size, 1 / (FS * oversample))
    in_band = (frequencies >= 300) & (frequencies <= 4_000)
    spectrum = np.zeros(frequencies.size, dtype=complex)
    count = np.count_nonzero(in_band)
    spectrum[in_band] = rng.normal(size=count) + 1j * rng.normal(size=count)
    source = np.fft.irfft(spectrum, grid.size)
    source *= 0.05 / np.std(source)
    truth = source[64 * oversample:(64 + length) * oversample:oversample]
    reference = truth + rng.normal(0, 0.0008, length)
    for ppm in (+100, -75):
        remote_positions = (index - 3.4) / (1 + ppm * 1e-6)
        remote = np.interp(remote_positions, grid, source) + rng.normal(0, 0.0008, length)
        yield ppm, reference, remote, truth


def _synthetic() -> dict:
    results = {}
    for ppm, reference, remote, truth in _synthetic_signals():
        result = _run_pair(reference, remote, truth)
        result["true_ppm"] = ppm
        aligned = [row for row in result["blocks"] if row["status"] == "aligned"]
        if aligned:
            errors = [row["predicted_lag_samples"] - (3.4 + ppm * 1e-6
                      * (row["start_sample"] + BLOCK // 2)) for row in aligned]
            result["aligned_center_lag_rmse_samples"] = float(np.sqrt(np.mean(np.square(errors))))
        results[f"{ppm:+d}_ppm"] = result
    return results


def _read_recording(raw: bytes) -> np.ndarray:
    with wave.open(io.BytesIO(raw), "rb") as wav:
        start = REAL_FIRST_BLOCK * BLOCK
        frames = REAL_BLOCKS * BLOCK
        if (wav.getframerate(), wav.getnchannels(), wav.getsampwidth()) != (FS, 1, 2):
            raise ValueError("expected a LibriWASN mono PCM16 recording at 16 kHz")
        if wav.getnframes() < start + frames:
            raise ValueError("LibriWASN recording is shorter than the experiment interval")
        wav.setpos(start)
        data = wav.readframes(frames)
    if len(data) != frames * 2:
        raise ValueError("truncated LibriWASN PCM data")
    return np.frombuffer(data, dtype="<i2").astype(float) / 32768


def _libriwasn_inputs(source: Path) -> tuple[dict, dict]:
    members = {device: device + "_ntlab_asn_ID_0.wav" for device in _DEVICES}
    raw = {}
    if source.is_dir():
        for device, name in members.items():
            with (source / name).open("rb") as stream:
                raw[device] = stream.read(_WAV_BYTES[device] + 1)
    else:
        with zipfile.ZipFile(source) as archive:
            for device, name in members.items():
                with archive.open(_WAV_MEMBER_ROOT + name) as stream:
                    raw[device] = stream.read(_WAV_BYTES[device] + 1)
    if any(len(raw[device]) != _WAV_BYTES[device] for device in members):
        raise ValueError("LibriWASN member sizes do not match the fixed session8 recordings")
    digests = {device: hashlib.sha256(data).hexdigest() for device, data in raw.items()}
    if digests != _WAV_SHA256:
        raise ValueError("LibriWASN member hashes do not match the fixed session8 recordings")
    samples = {device: _read_recording(data) for device, data in raw.items()}
    metadata = {
        "dataset_doi": "10.5281/zenodo.7960972",
        "archive": "LibriWASN_200_0L.zip",
        "archive_members": {device: _WAV_MEMBER_ROOT + name for device, name in members.items()},
        "wav_sha256": digests,
        "reference": _DEVICES[0],
        "interval_samples": [REAL_FIRST_BLOCK * BLOCK, (REAL_FIRST_BLOCK + REAL_BLOCKS) * BLOCK],
    }
    return samples, metadata


def _libriwasn(source: Path) -> dict:
    samples, metadata = _libriwasn_inputs(source)
    return {**metadata, "pairs": {device: _run_pair(samples[_DEVICES[0]], samples[device])
                                  for device in _DEVICES[1:]}}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--libriwasn", type=Path, help="official ZIP or extracted session8 record directory")
    args = parser.parse_args()
    result = {
        "sample_rate_hz": FS,
        "block_samples": BLOCK,
        "extra_buffer_samples": BUFFER,
        "dropped_blocks_zero_based": list(range(*GAP)),
        "synthetic": _synthetic(),
    }
    if args.libriwasn is not None:
        result["libriwasn"] = _libriwasn(args.libriwasn)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
