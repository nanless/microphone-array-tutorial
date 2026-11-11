"""Reproduce the course's LOCATA DICIT listening pair from extracted dev.zip members."""
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
from codes.chapters.ch05.core.beamforming import apply_beamformer, dsb_weights

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

    # Known source A and measured microphones give spherical-wave travel times.
    distances = np.linalg.norm(sources[0] - microphones, axis=1)
    delays = (distances - distances[REFERENCE]) / 343.0
    steering = np.exp(-2j * np.pi * frequencies[:, None] * delays)
    beam = apply_beamformer(spectra, dsb_weights(steering))
    after = istft(beam[None], n_fft=N_FFT, hop_length=HOP, length=waveforms.shape[1])[0]

    channels = np.stack([waveforms[REFERENCE], after])
    scaled = np.rint(channels * GAIN * 32768)
    if np.max(np.abs(scaled)) > 32767:
        raise ValueError("Common gain would clip")
    pcm = scaled.astype("<i2")
    assets = {name: wav_bytes(values) for name, values in zip(("before.wav", "after.wav"), pcm)}
    energy_ratio = {}
    for name, mask in (("A_only", solo_a), ("B_only", solo_b)):
        energy = np.sum(pcm[:, mask].astype(np.int64) ** 2, axis=1)
        energy_ratio[name] = {
            "samples": int(mask.sum()), "before_pcm_energy": int(energy[0]),
            "after_pcm_energy": int(energy[1]), "change_db": float(10 * np.log10(energy[1] / energy[0])),
        }
    manifest = {
        "record": "https://zenodo.org/records/3630471",
        "members_sha256": {RECORDING + name: digest for name, digest in MEMBERS.items()},
        "files_sha256": {name: sha256(data).hexdigest() for name, data in assets.items()},
        "sample_rate": FS, "samples": waveforms.shape[1], "reference_channel_1_based": REFERENCE + 1,
        "export_gain": GAIN, "fft": N_FFT, "hop": HOP, "sound_speed_m_s": 343.0,
        "sources_xyz_m_relative_to_array": sources.tolist(),
        "source_bearings_deg_from_positive_y": np.degrees(np.arctan2(sources[:, 0], sources[:, 1])).tolist(),
        "doa_scan": {"band_hz": [375, 3000], "grid_deg": angles.tolist(), "scores": scores.tolist(),
                     "peak_deg": int(angles[np.argmax(scores)]), "frames": int(frame_mask.sum())},
        "energy_ratio": energy_ratio,
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
