"""Compare audio-estimated steering and optical steering on LOCATA dev Task 2."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import numpy as np

from codes.chapters.ch02.core.spectral import istft, stft
from codes.chapters.ch03.core.geometry import plane_wave_steering
from codes.chapters.ch04.core.doa import srp_phat
from codes.chapters.ch05.core.beamforming import apply_beamformer, dsb_weights, mvdr_weights
from scripts.course_home import prepare_locata as locata


# LOCATA final-release documentation, Figure 6 and Table 1. The y coordinate
# is zero for every microphone in the DICIT local coordinate system.
MICROPHONES = np.column_stack((
    [0.96, 0.64, 0.32, 0.16, 0.08, 0.04, 0, 0.96, -0.04, -0.08, -0.16, -0.32, -0.64, -0.96, -0.96],
    np.zeros(locata.CHANNELS),
    [0, 0, 0, 0, 0, 0, 0, 0.32, 0, 0, 0, 0, 0, 0, 0.32],
))
A_SCAN = (124320, 135360)  # 2.59–2.82 s
SCORE_WINDOWS = {"A_only": (137280, 148800), "B_only": locata.WINDOWS["B_only"]}
SCAN_BAND = (375, 3000)
ANGLES_DEG = np.arange(-90, 91)
AZIMUTH_OFFSET_DEG = -20.0


def run(input_dir: Path) -> dict:
    """Read the fixed extracted recording and return powers relative to raw mic 7."""
    paths = {
        name: input_dir / locata.RECORDING / name for name in locata.MEMBERS
        if name != "benchmark2/position_source_loudspeaker4.txt"
    }
    for name, path in paths.items():
        if sha256(path.read_bytes()).hexdigest() != locata.MEMBERS[name]:
            raise ValueError(f"Unexpected LOCATA member: {name}")

    def read_table(name):
        return np.loadtxt(paths[name], skiprows=1)

    source = paths["dicit/audio_array_dicit.wav"]
    raw = subprocess.check_output(["ffmpeg", "-v", "error", "-i", str(source), "-f", "f64le", "-"])
    waveforms = np.frombuffer(raw, dtype="<f8").reshape(-1, locata.CHANNELS).T
    if waveforms.shape != (locata.CHANNELS, 302528) or not np.isfinite(waveforms).all():
        raise ValueError("Unexpected DICIT recording shape or values")

    activity_a, activity_b = [
        read_table(f"dicit/VAD_dicit_loudspeaker{index}.txt").astype(bool) for index in (3, 4)
    ]
    solo_a, solo_b = activity_a & ~activity_b, activity_b & ~activity_a
    for (start, stop), label in (
        (A_SCAN, solo_a), (SCORE_WINDOWS["A_only"], solo_a),
        (locata.CALIBRATION, solo_b), (SCORE_WINDOWS["B_only"], solo_b),
    ):
        if not label[start:stop].all():
            raise ValueError("Expected a single-speaker interval")

    input_spectra = stft(waveforms, n_fft=locata.N_FFT, hop_length=locata.HOP)
    frequencies = np.fft.rfftfreq(locata.N_FFT, 1 / locata.FS)
    centers = np.arange(input_spectra.shape[-1]) * locata.HOP
    scan_frames = locata.complete_frames(centers, A_SCAN)
    scan_bins = (frequencies >= SCAN_BAND[0]) & (frequencies <= SCAN_BAND[1])
    scan = srp_phat(input_spectra[:, scan_bins][:, :, scan_frames], frequencies[scan_bins],
                    MICROPHONES, np.radians(ANGLES_DEG))
    estimated_deg = int(ANGLES_DEG[np.argmax(scan)])

    # Fit one interference covariance and reuse it for every steering vector.
    covariance_frames = locata.complete_frames(centers, locata.CALIBRATION)
    training = input_spectra[:, :, covariance_frames]
    covariance = np.einsum("mft,nft->fmn", training, training.conj()) / training.shape[-1]
    mix = locata.raised_cosine(frequencies, 125, 250) * locata.raised_cosine(-frequencies, -8000, -6000)

    def plane(angle_deg):
        return plane_wave_steering(MICROPHONES, frequencies, np.radians(angle_deg), reference=locata.REFERENCE)

    def beam_outputs(steering):
        ds = dsb_weights(steering)
        mvdr = mvdr_weights(covariance, steering, relative_diagonal_loading=locata.LOADING)
        mvdr = ds + mix[:, None] * (mvdr - ds)
        return [
            istft(apply_beamformer(input_spectra, weights)[None], n_fft=locata.N_FFT,
                  hop_length=locata.HOP, length=waveforms.shape[-1])[0]
            for weights in (ds, mvdr)
        ]

    output = {"mic_7": waveforms[locata.REFERENCE]}

    def add_beams(name, steering):
        ds, mvdr = beam_outputs(steering)
        output[f"ds_{name}"] = ds
        output[f"mvdr_{name}"] = mvdr

    add_beams("estimated", plane(estimated_deg))
    add_beams("offset", plane(estimated_deg + AZIMUTH_OFFSET_DEG))

    # Only the reference directions below use motion capture.
    array = np.median(read_table("dicit/position_array_dicit.txt"), axis=0)
    target = np.median(read_table("benchmark2/position_source_loudspeaker3.txt"), axis=0)[6:9]
    local_target = (target - array[6:9]) @ array[12:21].reshape(3, 3)
    optical_deg = float(np.degrees(np.arctan2(local_target[0], local_target[1])))
    add_beams("optical_azimuth", plane(optical_deg))
    distances = np.linalg.norm(target - array[21:].reshape(locata.CHANNELS, 3), axis=1)
    delays = (distances - distances[locata.REFERENCE]) / 343.0
    add_beams("optical_position", np.exp(-2j * np.pi * frequencies[:, None] * delays))

    names = list(output)
    scaled = np.rint(np.stack(list(output.values())) * locata.GAIN * 32768)
    if not np.isfinite(scaled).all() or scaled.min() < -32768 or scaled.max() > 32767:
        raise ValueError("Cannot export the comparison to PCM16 at the common gain")
    pcm = scaled.astype("<i2")
    scored = {}
    for name, (start, stop) in SCORE_WINDOWS.items():
        band_power = locata.band_mean_square(pcm[:, start:stop], locata.EVALUATION_BAND)
        scored[name] = {
            "sample_range": [start, stop],
            "band_mean_square_fs2": dict(zip(names, band_power.tolist())),
            "band_db_vs_mic_7": dict(zip(names, (10 * np.log10(band_power / band_power[0])).tolist())),
        }
    return {
        "input_sha256": {locata.RECORDING + name: locata.MEMBERS[name] for name in paths},
        "scan_samples": list(A_SCAN), "scan_frames": int(scan_frames.sum()), "scan_band_hz": list(SCAN_BAND),
        "estimated_deg": estimated_deg, "optical_azimuth_deg": optical_deg,
        "offset_control_deg": estimated_deg + AZIMUTH_OFFSET_DEG,
        "interference_fit_samples": list(locata.CALIBRATION), "interference_fit_frames": int(covariance_frames.sum()),
        "fft": locata.N_FFT, "hop": locata.HOP, "sample_rate_hz": locata.FS,
        "relative_diagonal_loading": locata.LOADING, "pcm_export_gain": locata.GAIN,
        "score_band_hz": list(locata.EVALUATION_BAND), "score_windows": scored,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path, help="Directory containing the extracted dev/ tree")
    args = parser.parse_args()
    print(json.dumps(run(args.input_dir), indent=2, allow_nan=False))
