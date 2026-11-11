"""Reproduce the course's LOCATA DICIT comparison from extracted dev.zip members."""
import argparse
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import subprocess
import wave

import numpy as np

from codes.chapters.ch02.core.spectral import istft, stft
from codes.chapters.ch04.core.doa import srp_phat
from codes.chapters.ch05.core.beamforming import apply_beamformer, dsb_weights, mvdr_weights

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "assets/readme/cases/locata"
RECORDING = "dev/task2/recording1/"
MEMBERS = {
    "dicit/audio_array_dicit.wav": "76b949375cee3185983aaf7cd839840f6ba76a88f7bbea6d308d376f5813a9dc",
    "dicit/position_array_dicit.txt": "f5dd34f1534991596c343550eb8613f75eabcee6dd69ade7e46544055fe28cc9",
    "dicit/VAD_dicit_loudspeaker3.txt": "3aa6db7eed941930a08241cf2261669df4cd90db72236bb26867d2dfd6e61c09",
    "dicit/VAD_dicit_loudspeaker4.txt": "096d42819b9a6ceb5171fa43c5dd35ff5160f520de20774db3b318e73341797e",
    # The release carries the same global source tracks under each array.
    "benchmark2/position_source_loudspeaker3.txt": "c43f4b72d0fdad7c80ebaad866cf14a313820ec3648e9860f95dabca2c9bdfef",
    "benchmark2/position_source_loudspeaker4.txt": "68d502dbc1438ac3b0f0925fdedcb8170bcbc55cd7bdc3d36925127c5c2952c4",
}
FS = 48000
CHANNELS = 15
REFERENCE = 6
N_FFT = 1024
HOP = 512
GAIN = 4

# These sample ranges are half-open and are disjoint from the fit range.
CALIBRATION = (167040, 189120)  # 3.48–3.94 s, B only.
WINDOWS = {
    "A_only": (127200, 148800),  # 2.65–3.10 s
    "B_only": (195840, 224640),  # 4.08–4.68 s
    "overlap": (235200, 273600),  # 4.90–5.70 s
}
EVALUATION_BAND = (300, 3400)
LOADING = 0.1
FILES = ("before.wav", "after.wav", "mvdr.wav", "focus.wav")
POSTFILTER_STRENGTH = 0.5
POSTFILTER_GAIN_FLOOR = 0.25


def band_mean_square(pcm, band):
    """Hann-weighted mean square in an interior one-sided DFT band, in FS²."""
    length = pcm.shape[-1]
    window = np.hanning(length)
    freqs = np.fft.rfftfreq(length, 1 / FS)
    band_bins = (freqs >= band[0]) & (freqs <= band[1])
    bins = np.fft.rfft(pcm * (window / 32768), axis=-1)[..., band_bins]
    return 2 * np.sum(abs(bins)**2, axis=-1) / (length * (window @ window))


def raised_cosine(x, low, high):
    fraction = np.clip((x - low) / (high - low), 0, 1)
    return 0.5 - 0.5 * np.cos(np.pi * fraction)


def level_plot(pcm):
    """Time-aligned band levels for the three published PCM signals."""
    block, hop = 2400, 960
    starts = np.arange(0, pcm.shape[1] - block + 1, hop)
    frames = np.lib.stride_tricks.sliding_window_view(pcm, block, axis=1)[:, starts]
    power = band_mean_square(frames, EVALUATION_BAND)
    seconds = (starts + block / 2) / FS
    x0, x1, y0, y1 = 56, 602, 45, 209
    duration = pcm.shape[1] / FS
    floor, ceiling = -60, -15

    def x(t):
        return x0 + (x1 - x0) * t / duration

    def y(db):
        return y0 + (y1 - y0) * (ceiling - db) / (ceiling - floor)

    parts = [
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 620 254" '
        'role="img" aria-labelledby="level-title level-desc" '
        'font-family="system-ui,sans-serif" font-size="17">',
        '<title id="level-title">同一录音的三路输出电平</title>',
        '<desc id="level-desc">300 至 3400 赫兹的时序电平。浅蓝标出目标单独发声，'
        '浅橙标出干扰单独发声，浅紫标出重叠。干扰时段的蓝色 MVDR 加后滤波曲线低于对齐相加。</desc>',
        '<rect width="620" height="254" rx="12" fill="#f6f8fb"/>',
        '<text x="16" y="28" fill="#334a64">dBFS</text>',
    ]
    bands = (
        ("A_only", "#dce8fb", "A"),
        ("B_only", "#f6e4d5", "B"),
        ("overlap", "#e8e3f4", "A + B"),
    )
    for name, color, label in bands:
        start, stop = (point / FS for point in WINDOWS[name])
        parts.append(
            f'<rect x="{x(start):.2f}" y="{y0}" width="{x(stop)-x(start):.2f}" '
            f'height="{y1-y0}" fill="{color}"/>'
        )
        parts.append(
            f'<text x="{x((start+stop)/2):.2f}" y="28" text-anchor="middle" '
            f'fill="#334a64" font-weight="700">{label}</text>'
        )
    for db in (-60, -45, -30, -15):
        parts.append(
            f'<path d="M{x0} {y(db):.2f}H{x1}" stroke="#c2cfdb"/>'
            f'<text x="{x0-7}" y="{y(db)+5:.2f}" text-anchor="end" fill="#334a64">{db}</text>'
        )
    for second in range(7):
        parts.append(
            f'<path d="M{x(second):.2f} {y1}v5" stroke="#778ba0"/>'
            f'<text x="{x(second):.2f}" y="232" text-anchor="middle" fill="#334a64">{second}</text>'
        )
    for levels, color, dash in zip(power[[0, 1, 3]], ("#586b80", "#aa5d28", "#234b91"), ("3 3", "7 4", "")):
        db = 10 * np.log10(np.maximum(levels, 10 ** (floor / 10)))
        points = " ".join(f"{x(t):.2f},{y(d):.2f}" for t, d in zip(seconds, np.minimum(db, ceiling)))
        parts.append(
            f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2.1" '
            f'stroke-dasharray="{dash}" stroke-linejoin="round"/>'
        )
    parts.extend(['<text x="607" y="232" text-anchor="end" fill="#334a64">秒</text>', "</svg>"])
    return ("\n".join(parts) + "\n").encode()


def render_beams(spectra, frequencies, steering, other, centers, solo_b, length):
    """Fit fixed MVDR weights only on complete B-only calibration frames."""
    start, stop = CALIBRATION
    if not np.all(solo_b[start:stop]):
        raise ValueError("LOCATA calibration is no longer B-only")
    calibration = complete_frames(centers, CALIBRATION)
    training = spectra[:, :, calibration]
    covariance = np.einsum("mft,nft->fmn", training, training.conj()) / training.shape[-1]

    ds = dsb_weights(steering)
    # Preserve DS at the edges for the Chapter 5 reference comparison.
    mix = raised_cosine(frequencies, 125, 250) * raised_cosine(-frequencies, -8000, -6000)
    optimal = mvdr_weights(covariance, steering, relative_diagonal_loading=LOADING)
    mvdr = ds + mix[:, None] * (optimal - ds)
    weight_response = np.sum(mvdr.conj() * steering, axis=1)
    if not np.allclose(weight_response, 1, rtol=0, atol=1e-12):
        raise ValueError("MVDR and DS blend must preserve the design response")

    # Project the B beam off the A steering vector. It is a B monitor with
    # a modelled zero response to A, including where A and B are collinear.
    monitor = dsb_weights(other)
    monitor -= ds * np.sum(steering.conj() * monitor, axis=1, keepdims=True)
    a_residual = np.max(abs(np.sum(monitor.conj() * steering, axis=1)))
    if a_residual > 1e-12:
        raise ValueError("Interference monitor must block the design target")
    target = apply_beamformer(spectra, optimal)
    interference = apply_beamformer(spectra, monitor)
    target_power = local_power(target)
    interference_power = local_power(interference)
    ratio = target_power / np.maximum(target_power + POSTFILTER_STRENGTH * interference_power, 1e-40)
    focused = target * np.maximum(POSTFILTER_GAIN_FLOOR, ratio)
    beams = (apply_beamformer(spectra, ds), apply_beamformer(spectra, mvdr), focused)
    outputs = [
        istft(beam[None], n_fft=N_FFT, hop_length=HOP, length=length)[0] for beam in beams
    ]
    diagnostics = {
        "calibration_samples": list(CALIBRATION),
        "calibration_frames": int(calibration.sum()),
        "loading_relative_to_mean_diagonal": LOADING,
        "mvdr_full_weight_band_hz": [250, 6000],
        "mvdr_taper_outer_edges_hz": [125, 8000],
        "min_model_white_noise_gain_db": float(-10 * np.log10(np.max(np.sum(abs(mvdr)**2, axis=1)))),
        "max_model_response_error": float(np.max(abs(weight_response - 1))),
        "strong_demo": {
            "model_target_blocking_error": float(a_residual),
            "postfilter_strength": POSTFILTER_STRENGTH,
            "postfilter_gain_floor": POSTFILTER_GAIN_FLOOR,
            "power_smoothing": [0.25, 0.5, 0.25],
            "power_smoothing_edge": "zero padding",
            "look_ahead_frames": 1,
        },
    }
    return outputs, diagnostics


def complete_frames(centers, interval):
    start, stop = interval
    return (centers - N_FFT // 2 >= start) & (centers + N_FFT // 2 <= stop)


def local_power(spectrum):
    power = np.pad(abs(spectrum)**2, ((0, 0), (1, 1)))
    return 0.25 * power[:, :-2] + 0.5 * power[:, 1:-1] + 0.25 * power[:, 2:]


def wav_bytes(pcm):
    stream = BytesIO()
    with wave.open(stream, "wb") as writer:
        writer.setparams((1, 2, FS, len(pcm), "NONE", "not compressed"))
        writer.writeframes(pcm.tobytes())
    return stream.getvalue()


def generate(input_dir):
    paths = {name: input_dir / RECORDING / name for name in MEMBERS}
    for name, path in paths.items():
        if sha256(path.read_bytes()).hexdigest() != MEMBERS[name]:
            raise ValueError(f"Unexpected LOCATA member: {name}")

    def read_table(name):
        return np.loadtxt(paths[name], skiprows=1)

    def median_row(name):
        return np.median(read_table(name), axis=0)

    source_file = paths["dicit/audio_array_dicit.wav"]
    raw = subprocess.check_output(
        ["ffmpeg", "-v", "error", "-i", str(source_file), "-f", "f64le", "-"]
    )
    waveforms = np.frombuffer(raw, dtype="<f8").reshape(-1, CHANNELS).T
    if waveforms.shape != (CHANNELS, 302528) or not np.all(np.isfinite(waveforms)):
        raise ValueError("Unexpected DICIT recording shape or values")
    geometry = median_row("dicit/position_array_dicit.txt")
    center = geometry[6:9]
    microphones = geometry[21:].reshape(CHANNELS, 3) - center
    sources = np.stack([median_row(f"benchmark2/position_source_loudspeaker{i}.txt")[6:9] - center
                        for i in (3, 4)])
    activity = [read_table(f"dicit/VAD_dicit_loudspeaker{i}.txt").astype(bool) for i in (3, 4)]
    solo_a = activity[0] & ~activity[1]
    solo_b = activity[1] & ~activity[0]
    if len(solo_a) != waveforms.shape[1] or min(solo_a.sum(), solo_b.sum()) == 0:
        raise ValueError("Expected two non-empty, sample-aligned activity masks")

    spectra = stft(waveforms, n_fft=N_FFT, hop_length=HOP)
    frequencies = np.fft.rfftfreq(N_FFT, 1 / FS)
    angles = np.arange(-90, 91)
    frequency_mask = (frequencies >= 375) & (frequencies <= 3000)
    centers = np.arange(spectra.shape[-1]) * HOP
    frame_mask = (centers < len(solo_a)) & solo_a[np.minimum(centers, len(solo_a) - 1)]
    scores = srp_phat(spectra[:, frequency_mask][:, :, frame_mask], frequencies[frequency_mask],
                      microphones, np.radians(angles))

    localization = {}
    for name in ("A_only", "B_only"):
        frames = complete_frames(centers, WINDOWS[name])
        selected = spectra[:, frequency_mask][:, :, frames]
        curve = srp_phat(selected, frequencies[frequency_mask], microphones, np.radians(angles))
        localization[name] = {
            "sample_range": list(WINDOWS[name]), "frames": int(frames.sum()),
            "peak_deg": int(angles[np.argmax(curve)]), "scores": curve.tolist(),
        }

    # Measured source and microphone coordinates give spherical-wave travel times.
    distances = np.linalg.norm(sources[:, None, :] - microphones[None, :, :], axis=2)
    delays = (distances - distances[:, [REFERENCE]]) / 343.0
    steering = np.exp(-2j * np.pi * frequencies[:, None, None] * delays)
    outputs, diagnostics = render_beams(spectra, frequencies, steering[:, 0], steering[:, 1],
                                        centers, solo_b, waveforms.shape[1])
    scaled = np.rint(np.stack([waveforms[REFERENCE], *outputs]) * GAIN * 32768)
    if np.max(np.abs(scaled)) > 32767:
        raise ValueError("Common gain would clip")
    pcm = scaled.astype("<i2")
    assets = {name: wav_bytes(values) for name, values in zip(FILES, pcm)}
    assets["levels.svg"] = level_plot(pcm)

    energy_ratio = {}
    for name, mask in (("A_only", solo_a), ("B_only", solo_b)):
        # Preserve the earlier raw-versus-DS full-mask results.
        energy = np.sum(pcm[:2, mask].astype(np.int64) ** 2, axis=1)
        energy_ratio[name] = {
            "samples": int(mask.sum()), "before_pcm_energy": int(energy[0]),
            "after_pcm_energy": int(energy[1]), "change_db": float(10 * np.log10(energy[1] / energy[0])),
        }
    comparisons = {}
    labels = {"A_only": solo_a, "B_only": solo_b, "overlap": activity[0] & activity[1]}
    for name, (start, stop) in WINDOWS.items():
        if not (stop <= CALIBRATION[0] or start >= CALIBRATION[1]) or not labels[name][start:stop].all():
            raise ValueError(f"Invalid disjoint LOCATA comparison window: {name}")
        region = pcm[:, start:stop]
        energy = np.sum(region.astype(np.int64)**2, axis=-1)
        band_power = band_mean_square(region, EVALUATION_BAND)
        comparisons[name] = {
            "sample_range": [start, stop],
            "pcm_sum_squares": dict(zip(FILES, energy.tolist())),
            "full_band_db_vs_reference": dict(zip(FILES, (10 * np.log10(energy / energy[0])).tolist())),
            "band_mean_square_fs2": dict(zip(FILES, band_power.tolist())),
            "band_db_vs_reference": dict(zip(FILES, (10 * np.log10(band_power / band_power[0])).tolist())),
        }
    implementation_sources = (
        Path(__file__).resolve(),
        ROOT / "codes/chapters/ch02/core/conventions.py",
        ROOT / "codes/chapters/ch02/core/spectral.py",
        ROOT / "codes/chapters/ch03/core/geometry.py",
        ROOT / "codes/chapters/ch04/core/covariance.py",
        ROOT / "codes/chapters/ch04/core/doa.py",
        ROOT / "codes/chapters/ch05/core/beamforming.py",
    )
    manifest = {
        "record": "https://zenodo.org/records/3630471",
        "source_sha256": {
            path.relative_to(ROOT).as_posix(): sha256(path.read_bytes()).hexdigest()
            for path in implementation_sources
        },
        "members_sha256": {RECORDING + name: digest for name, digest in MEMBERS.items()},
        "files_sha256": {name: sha256(data).hexdigest() for name, data in assets.items()},
        "sample_rate": FS, "samples": waveforms.shape[1], "reference_channel_1_based": REFERENCE + 1,
        "export_gain": GAIN, "fft": N_FFT, "hop": HOP, "sound_speed_m_s": 343.0,
        "sources_xyz_m_relative_to_array": sources.tolist(),
        "source_bearings_deg_from_positive_y": np.degrees(np.arctan2(sources[:, 0], sources[:, 1])).tolist(),
        "doa_scan": {"band_hz": [375, 3000], "grid_deg": angles.tolist(), "scores": scores.tolist(),
                     "peak_deg": int(angles[np.argmax(scores)]), "frames": int(frame_mask.sum())},
        "localization_windows": localization,
        "energy_ratio": energy_ratio,
        "mvdr": diagnostics,
        "comparison": {"band_hz": list(EVALUATION_BAND), "windows": comparisons},
        "level_plot": {"window_samples": 2400, "hop_samples": 960, "floor_dbfs": -60,
                       "ceiling_dbfs": -15},
    }
    assets["MANIFEST.json"] = (json.dumps(manifest, indent=2, allow_nan=False) + "\n").encode()
    return assets


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path, help="Directory containing the extracted dev/ tree")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    for name, content in generate(args.input_dir).items():
        path = OUT / name
        if args.check:
            if path.read_bytes() != content:
                raise SystemExit(f"Published LOCATA asset differs: {name}")
        else:
            path.write_bytes(content)


if __name__ == "__main__":
    main()
