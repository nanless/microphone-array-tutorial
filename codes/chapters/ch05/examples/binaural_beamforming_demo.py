"""Known-component controls for two-output beamforming; no STFT or device model."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from codes.chapters.ch00.core.audio_samples import pcm16_bytes, read_pcm16
from codes.chapters.ch00.io_contracts import (
    same_metadata, strict_json_loads, validate_asset_directory, validate_parent_chain,
)
from codes.chapters.ch05.core.beamforming import lcmv_weights, mvdr_weights

ROOT = Path(__file__).resolve().parents[4]
OUTPUT = ROOT / "codes/chapters/ch05/binaural_beamforming_audio"
ROLES = ("input", "bmvdr", "blcmv", "null")
MEMBERS = {f"binaural_{role}.wav" for role in ROLES} | {"MANIFEST.json"}
SOURCES = (
    "codes/chapters/ch05/examples/binaural_beamforming_demo.py",
    "codes/chapters/ch05/core/beamforming.py",
    "codes/chapters/ch04/core/covariance.py",
    "codes/chapters/ch02/core/conventions.py",
    "codes/chapters/ch00/core/audio_samples.py",
    "codes/chapters/ch00/io_contracts.py",
)
FS, LENGTH, FADE = 16000, 32000, 320
START, STOP = 2400, 29600
FREQUENCIES = (500, 1500, 750, 1000, 1250)
AMPLITUDES = (0.12, 0.12, 0.08, 0.08, 0.08)


def fixture_weights():
    a = np.array([1., 1., 1.])
    b = np.array([1., -1., 1.])
    r = np.eye(3)
    reference = r[:, :2]
    one = mvdr_weights(r, a)
    c = np.column_stack((a, b))

    def constrained(eta):
        return np.column_stack([
            lcmv_weights(r, c, np.conj([a[i], eta * b[i]])) for i in range(2)
        ])

    weights = {
        "input": reference, "bmvdr": np.column_stack((one, one)),
        "blcmv": constrained(1 / 3), "null": constrained(0),
    }
    return a, b, weights


def _responses(a, b, weights):
    return {
        role: {
            "weights": w.real.T.tolist(),
            "target_response": (w.conj().T @ a).real.tolist(),
            "interferer_response": (w.conj().T @ b).real.tolist(),
            "unit_background_output_covariance": (w.conj().T @ w).real.tolist(),
        }
        for role, w in weights.items()
    }


def run_exercises():
    a, b, weights = fixture_weights()
    report = _responses(a, b, weights)
    # Nonreal directional constraint makes accidental conjugation observable.
    b_complex = np.array([1, 1j, 1])
    c = np.column_stack((a, b_complex))
    right = lcmv_weights(np.eye(3), c, np.conj([1, 1j / 3]))
    return {
        "E05-25": {role: row["interferer_response"] for role, row in report.items()},
        "E05-26": {role: row["unit_background_output_covariance"] for role, row in report.items()},
        "E05-27": {"null_weights": report["null"]["weights"],
                   "complex_right_response_real_imag": [
                       float(np.vdot(right, b_complex).real),
                       float(np.vdot(right, b_complex).imag),
                   ]},
    }


def _pcm_scores(samples):
    pcm = np.rint(samples[START:STOP] * 32768).astype(np.int64)
    energy = [int(pcm[:, ear] @ pcm[:, ear]) for ear in range(2)]
    d = (STOP - START) * 32768**2
    return {"integer_energy_each_ear": energy,
            "integer_cross_sum": int(pcm[:, 0] @ pcm[:, 1]),
            "integer_mean_square_denominator_per_ear": d}


def _tone_scores(x, basis):
    window = x[START:STOP]
    # The five tones are orthogonal on this fixed integer-cycle interval.
    return (2 / len(window) * (basis[:, START:STOP] @ window)).tolist()


def prepare_assets():
    a, b, weights = fixture_weights()
    # Verify every real weight against the independent rational solutions before export.
    expected = {
        "input": [[1, 0, 0], [0, 1, 0]],
        "bmvdr": [[1/3, 1/3, 1/3], [1/3, 1/3, 1/3]],
        "blcmv": [[1/3, 1/3, 1/3], [1/6, 2/3, 1/6]],
        "null": [[1/4, 1/2, 1/4], [1/4, 1/2, 1/4]],
    }
    n = np.arange(LENGTH)
    basis = np.cos((2 * np.pi / FS) * np.outer(FREQUENCIES, n))
    ramp = np.sin(np.linspace(0, np.pi / 2, FADE))**2
    envelope = np.ones(LENGTH)
    envelope[:FADE], envelope[-FADE:] = ramp, ramp[::-1]
    transfer = np.column_stack((a, b, np.eye(3)))
    observations = transfer @ (np.asarray(AMPLITUDES)[:, None] * basis)
    report = _responses(a, b, weights)
    payloads = {}
    for role in ROLES:
        w = weights[role]
        if not np.allclose(w.T, expected[role], rtol=0, atol=1e-14):
            raise ValueError(f"unexpected binaural control weights: {role}")
        x = (w.conj().T @ observations).real.T * envelope[:, None]
        filename = f"binaural_{role}.wav"
        payload = pcm16_bytes(x.T, FS)
        rate, pcm = read_pcm16(payload)
        if rate != FS or pcm.shape != (2, LENGTH):
            raise ValueError("unexpected stereo PCM format")
        pcm = pcm.T
        error = float(np.max(np.abs(x - pcm)))
        if error > 0.5 / 32768 + 1e-15:
            raise ValueError("binaural control exceeded PCM16 rounding bound")
        payloads[filename] = payload
        report[role].update(
            file=filename, sha256=sha256(payload).hexdigest(),
            peak_quantization_error=error,
            float_signed_cosine_amplitude_by_tone=_tone_scores(x, basis),
            pcm_signed_cosine_amplitude_by_tone=_tone_scores(pcm, basis),
            actual_pcm_total=_pcm_scores(pcm),
        )
    return payloads, {
        "schema_version": 1, "exercises": ["E05-25", "E05-26", "E05-27"],
        "origin": "deterministic five-tone component control, not measured acoustics",
        "sample_rate_hz": FS, "samples_per_channel": LENGTH,
        "channel_order": ["left", "right"], "common_export_gain": 1.0,
        "target_transfer": a.tolist(), "interferer_transfer": b.tolist(),
        "frequencies_hz": list(FREQUENCIES), "amplitudes": list(AMPLITUDES),
        "tone_order": ["target", "interferer", "background_mic1", "background_mic2", "background_mic3"],
        "window_samples": [START, STOP], "fade_samples_each_end": FADE,
        "sample_format": "signed little endian PCM16, ties to even, decode scale 32768",
        "source_sha256": {
            path: sha256(validate_parent_chain(ROOT / path).read_bytes()).hexdigest() for path in SOURCES
        },
        "measurements": report,
    }


def check_assets(output=OUTPUT):
    output = validate_asset_directory(output, MEMBERS, check=True)
    payloads, expected = prepare_assets()
    observed = strict_json_loads((output / "MANIFEST.json").read_bytes())
    if not same_metadata(observed, expected):
        raise ValueError("binaural control manifest differs from current source and replay")
    for filename, payload in payloads.items():
        if (output / filename).read_bytes() != payload:
            raise ValueError(f"binaural control PCM differs from replay: {filename}")
    return observed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--write", action="store_true")
    action.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        print(f"Verified {len(check_assets()['measurements'])} stereo WAVs")
    elif args.write:
        output = validate_asset_directory(OUTPUT, MEMBERS, check=False)
        payloads, manifest = prepare_assets()
        payloads["MANIFEST.json"] = (json.dumps(manifest, indent=2, allow_nan=False) + "\n").encode()
        output.mkdir(exist_ok=True)
        for filename, payload in payloads.items():
            (output / filename).write_bytes(payload)
        check_assets(output)
    else:
        print(json.dumps(run_exercises(), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
