"""Generate six synthetic wearable-array WAVs; --check replays strictly without writing."""
from __future__ import annotations

if __name__ == "__main__" and not __package__:
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[4]))

import argparse
import hashlib
import io
import json
from pathlib import Path
import platform
import wave

import numpy as np

from codes.chapters.ch00.core.audio_samples import pcm16_bytes, read_pcm16
from codes.chapters.ch00.io_contracts import (
    same_metadata,
    strict_json_loads,
    validate_asset_directory,
    validate_parent_chain,
)
from codes.chapters.ch14.core.wearable_audio import (
    CHANNEL_COUNTS,
    FILE_NAMES,
    LIMITS,
    SAMPLE_RATE,
    SAMPLES,
    measure_signal,
    parameters,
    run_experiment,
)

ROOT = Path(__file__).resolve().parents[4]
DEFAULT_OUTPUT = ROOT / "codes/chapters/ch14/wearable_audio"
SOURCE_PATHS = (
    "codes/chapters/ch14/core/wearable_fusion.py",
    "codes/chapters/ch14/core/wearable_audio.py",
    "codes/chapters/ch14/examples/generate_wearable_audio.py",
    "codes/chapters/ch02/core/conventions.py",
    "codes/chapters/ch00/core/audio_samples.py",
    "codes/chapters/ch00/io_contracts.py",
)

_FIXED_FILES = {
    name: name + ".wav"
    for name in (
        "wearable_clean_ref",
        "wearable_air_wind_2ch",
        "wearable_inear_occluded",
        "wearable_inear_equalized",
        "wearable_fused_output",
        "wearable_hearthrough_comb",
    )
}
_FIXED_CHANNELS = {
    "wearable_clean_ref": 1,
    "wearable_air_wind_2ch": 2,
    "wearable_inear_occluded": 1,
    "wearable_inear_equalized": 1,
    "wearable_fused_output": 1,
    "wearable_hearthrough_comb": 1,
}
_FIXED_LIMITS = (
    "Deterministic dual-tone target (250 Hz and 4000 Hz) with finite-window orthogonal "
    "disturbance tones, known occlusion gains (2.0x at 250 Hz, 0.1x at 4000 Hz), and "
    "exact 32-sample (2 ms) causal delay; not speech, real wind turbulence, measured "
    "ear-canal transfer functions, blind RTF estimation, or a subjective listening study. "
    "No per-file normalization, clipping, dither, or post-hoc gain/delay fitting."
)


def _fixed_parameters() -> dict:
    return {
        "exercise_ids": ["E14-12", "E14-13", "E14-14", "E14-16"],
        "sample_rate_hz": 16000,
        "samples_per_channel": 32000,
        "duration_s": 2.0,
        "fade_samples_each_end": 1600,
        "envelope": "0.5*(1-cos(pi*r/1600)) on [0,1600) and [30400,32000), 1.0 on [1600,30400)",
        "scoring_interval_samples": [1600, 30400],
        "scoring_samples_per_channel": 28800,
        "target_frequencies_hz": [250, 4000],
        "target_amplitudes": [0.10, 0.10],
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
            "low_air_inear": [1.0 / 401.0, 200.0 / 401.0],
            "high_air_inear": [80.0 / 81.0, 10.0 / 81.0],
        },
        "hearthrough_delay_samples": 32,
        "hearthrough_delay_ms": 2.0,
        "hearthrough_gains": [1.0, 1.0],
        "common_export_gain": 1.0,
    }


def _numeric_metadata(actual, expected) -> bool:
    """Exact structures/types; allow tiny roundoff only in independently computed floats."""
    if type(actual) is not type(expected):
        return False
    if isinstance(expected, dict):
        return actual.keys() == expected.keys() and all(
            _numeric_metadata(actual[k], v) for k, v in expected.items()
        )
    if isinstance(expected, list):
        return len(actual) == len(expected) and all(
            _numeric_metadata(a, b) for a, b in zip(actual, expected)
        )
    if type(expected) is float:
        return np.isfinite(actual) and bool(
            np.isclose(actual, expected, rtol=2e-12, atol=1e-13)
        )
    return actual == expected


def _validate_fixed_contract() -> None:
    if (
        type(SAMPLE_RATE) is not int
        or SAMPLE_RATE != 16000
        or type(SAMPLES) is not int
        or SAMPLES != 32000
        or not same_metadata(FILE_NAMES, _FIXED_FILES)
        or not same_metadata(CHANNEL_COUNTS, _FIXED_CHANNELS)
        or not same_metadata(LIMITS, _FIXED_LIMITS)
        or not _numeric_metadata(parameters(), _fixed_parameters())
    ):
        raise ValueError("wearable audio fixed parameters/types/limits differ")


def _fixed_waveforms() -> tuple[dict[str, np.ndarray], dict[str, dict[str, np.ndarray]]]:
    """Independently reconstruct the six declared waveforms and components."""
    n = np.arange(32000, dtype=np.float64)
    env = np.ones(32000, dtype=np.float64)
    env[:1600] = 0.5 * (1.0 - np.cos(np.pi * np.arange(1600) / 1600.0))
    env[30400:] = 0.5 * (1.0 - np.cos(np.pi * (31999 - np.arange(30400, 32000)) / 1600.0))

    s_low = env * (0.10 * np.cos(2.0 * np.pi * 250.0 * n / 16000.0))
    s_high = env * (0.10 * np.cos(2.0 * np.pi * 4000.0 * n / 16000.0))
    s_clean = s_low + s_high

    w0 = env * (0.20 * np.cos(2.0 * np.pi * 125.0 * n / 16000.0))
    w1 = env * (0.20 * np.cos(2.0 * np.pi * 375.0 * n / 16000.0))
    v0 = env * (
        0.02 * np.cos(2.0 * np.pi * 3500.0 * n / 16000.0)
        + 0.01 * np.sin(2.0 * np.pi * 3500.0 * n / 16000.0)
    )
    v1 = env * (
        0.02 * np.cos(2.0 * np.pi * 4500.0 * n / 16000.0)
        + 0.01 * np.sin(2.0 * np.pi * 4500.0 * n / 16000.0)
    )
    air = np.vstack([s_clean + w0 + v0, s_clean + w1 + v1])

    vin_low = env * (0.02 * np.cos(2.0 * np.pi * 175.0 * n / 16000.0))
    vin_high = env * (0.02 * np.cos(2.0 * np.pi * 5000.0 * n / 16000.0))
    in_low = 2.0 * s_low + vin_low
    in_high = 0.1 * s_high + vin_high
    inear = in_low + in_high
    eq = 0.5 * in_low + 10.0 * in_high

    fused = (1.0 / 401.0) * (s_low + w0) + (200.0 / 401.0) * in_low + (80.0 / 81.0) * (
        s_high + v0
    ) + (10.0 / 81.0) * in_high
    fused_clean = s_low + s_high
    fused_noise = (
        (1.0 / 401.0) * w0
        + (200.0 / 401.0) * vin_low
        + (80.0 / 81.0) * v0
        + (10.0 / 81.0) * vin_high
    )

    delayed = np.zeros_like(s_clean)
    delayed[32:] = s_clean[:-32]
    ht = s_clean + delayed

    waves = {
        "wearable_clean_ref": s_clean[None, :],
        "wearable_air_wind_2ch": air,
        "wearable_inear_occluded": inear[None, :],
        "wearable_inear_equalized": eq[None, :],
        "wearable_fused_output": fused[None, :],
        "wearable_hearthrough_comb": ht[None, :],
    }
    comps = {
        "wearable_clean_ref": {"clean": s_clean, "noise": np.zeros_like(s_clean)},
        "wearable_air_wind_2ch": {"clean": s_clean, "noise": w0 + v0},
        "wearable_inear_occluded": {"clean": 2.0 * s_low + 0.1 * s_high, "noise": vin_low + vin_high},
        "wearable_inear_equalized": {"clean": s_clean, "noise": 0.5 * vin_low + 10.0 * vin_high},
        "wearable_fused_output": {"clean": fused_clean, "noise": fused_noise},
        "wearable_hearthrough_comb": {"clean": ht, "noise": np.zeros_like(s_clean)},
    }
    return waves, comps


def _independent_score(x: np.ndarray, ref: np.ndarray, comp: dict | None = None, *, pcm: bool = False) -> dict:
    win = x[:, 1600:30400]
    ref_win = ref[0, 1600:30400]
    idx = np.arange(1600, 30400, dtype=np.float64)
    proj = [
        2.0 / 28800.0 * np.sum(win[0] * np.exp(-2j * np.pi * f * idx / 16000.0))
        for f in (250.0, 4000.0)
    ]
    mse = float(np.mean((win[0] - ref_win) ** 2))
    ref_pow = float(np.mean(ref_win**2))
    out = {
        "scoring_samples_per_channel": 28800,
        "peak_all_samples": float(np.max(np.abs(x))),
        "mean_square_by_channel": [float(np.mean(c**2)) for c in win],
        "tone_projections_ch0_real_imag": [[float(z.real), float(z.imag)] for z in proj],
        "reference_mse_ch0": mse,
        "reference_nmse_ch0": float(mse / ref_pow),
    }
    if not pcm and comp is not None:
        cw = comp["clean"][1600:30400]
        nw = comp["noise"][1600:30400]
        dw = cw - ref_win
        out["float_decomposition_ch0"] = {
            "target_distortion_power": float(np.mean(dw**2)),
            "residual_noise_power": float(np.mean(nw**2)),
            "cross_term": float(2.0 * np.mean(dw * nw)),
            "reconstruction_max_abs_error": float(np.max(np.abs(win[0] - cw - nw))),
        }
    if pcm:
        xi = np.rint(x * 32768.0).astype(np.int64)[:, 1600:30400]
        ri = np.rint(ref[0] * 32768.0).astype(np.int64)[1600:30400]
        sq = [sum(int(v) * int(v) for v in ch) for ch in xi]
        err_sq = sum(int(a - b) * int(a - b) for a, b in zip(xi[0], ri))
        ref_sq = sum(int(v) * int(v) for v in ri)
        out["pcm_integer_metrics"] = {
            "pcm_decode_divisor": 32768,
            "integer_power_sample_denominator": 28800,
            "integer_power_scale_denominator": 32768**2,
            "integer_scored_squared_sum_by_channel": sq,
            "integer_reference_error_squared_sum_ch0": err_sq,
            "integer_reference_squared_sum": ref_sq,
            "integer_nmse_ch0": float(err_sq / ref_sq),
        }
    return out


def _validate_model(report: dict, arrays: dict[str, np.ndarray]) -> None:
    fixed_waves, fixed_comps = _fixed_waveforms()
    if (
        type(report) is not dict
        or not _numeric_metadata(report.get("parameters"), _fixed_parameters())
        or not same_metadata(report.get("limits"), _FIXED_LIMITS)
    ):
        raise ValueError("wearable audio model fixed metadata/types/limits differ")
    if type(arrays) is not dict or arrays.keys() != fixed_waves.keys():
        raise ValueError("wearable audio fixture members differ")
    ref = fixed_waves["wearable_clean_ref"]
    for key, expected in fixed_waves.items():
        arr = arrays[key]
        if (
            type(arr) is not np.ndarray
            or arr.dtype != np.dtype("float64")
            or arr.shape != expected.shape
            or not np.all(np.isfinite(arr))
            or not np.allclose(arr, expected, rtol=1e-14, atol=1e-15)
        ):
            raise ValueError("wearable audio fixed waveform/type/shape differs: " + key)
        indep = _independent_score(arr, ref, fixed_comps[key], pcm=False)
        if not _numeric_metadata(report["float_measurements"][key], indep):
            raise ValueError("wearable audio independent float score differs: " + key)


def _sha(blob: bytes) -> str:
    return hashlib.sha256(blob).hexdigest()


def source_digests() -> dict[str, str]:
    return {path: _sha(validate_parent_chain(ROOT / path).read_bytes()) for path in SOURCE_PATHS}


def prepare_assets() -> tuple[dict, dict[str, bytes]]:
    _validate_fixed_contract()
    report, arrays = run_experiment()
    _validate_model(report, arrays)

    contents: dict[str, bytes] = {}
    files: dict[str, dict] = {}
    decoded: dict[str, np.ndarray] = {}

    for key, arr in arrays.items():
        blob = pcm16_bytes(arr, SAMPLE_RATE)
        rate, pcm = read_pcm16(blob)
        channels = _FIXED_CHANNELS[key]
        with wave.open(io.BytesIO(blob)) as reader:
            if (
                reader.getframerate(),
                reader.getnchannels(),
                reader.getnframes(),
                reader.getsampwidth(),
                reader.getcomptype(),
            ) != (16000, channels, 32000, 2, "NONE"):
                raise ValueError("wearable audio in-memory WAV format differs: " + key)
            integers = np.frombuffer(reader.readframes(32000), dtype="<i2").reshape(-1, channels).T
        if (
            type(rate) is not int
            or rate != 16000
            or not np.array_equal(integers, np.rint(arr * 32768.0))
            or type(pcm) is not np.ndarray
            or pcm.dtype != np.dtype("float64")
            or not np.array_equal(pcm, integers.astype(np.float64) / 32768.0)
        ):
            raise ValueError("wearable audio independent PCM codec check failed: " + key)
        filename = FILE_NAMES[key]
        contents[filename] = blob
        decoded[key] = pcm
        files[filename] = {
            "channels": channels,
            "samples_per_channel": SAMPLES,
            "sample_rate_hz": rate,
            "sha256": _sha(blob),
        }

    pcm_ref = decoded["wearable_clean_ref"]
    samples: dict[str, dict] = {}
    for key, arr in arrays.items():
        filename = FILE_NAMES[key]
        pcm = decoded[key]
        pcm_scores = measure_signal(pcm, key, reference=pcm_ref, pcm=True)
        indep_pcm = _independent_score(pcm, pcm_ref, pcm=True)
        if not _numeric_metadata(pcm_scores, indep_pcm):
            raise ValueError("wearable audio independent PCM score differs: " + key)
        samples[key] = {
            "file": filename,
            "analytic_expectations": report["analytic_expectations"][key],
            "float_measurements": report["float_measurements"][key],
            "pcm_measurements": pcm_scores,
            "quantization_max_abs_error": float(np.max(np.abs(arr - pcm))),
        }

    manifest = {
        "schema_version": 1,
        "origin": "original deterministic dual-tone wearable array and cross-modal fusion controls",
        "sample_rate_hz": SAMPLE_RATE,
        "samples_per_channel": SAMPLES,
        "common_export_gain": 1.0,
        "source_sha256": source_digests(),
        "parameters": report["parameters"],
        "files": files,
        "samples": samples,
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "system": platform.system(),
            "machine": platform.machine(),
        },
        "pcm_encode_multiplier": 32768,
        "pcm_decode_divisor": 32768,
        "quantization_max_abs_error_bound": 0.5 / 32768.0,
        "pcm": "signed little-endian PCM16, nearest-even rounding; no dither, saturation or per-file gain",
        "limits": LIMITS,
        "listening": "Start at low device volume; no autoplay or subjective listening study.",
    }
    contents["MANIFEST.json"] = (
        json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    ).encode("utf-8")
    return manifest, contents


def check_assets(directory: Path | str = DEFAULT_OUTPUT) -> dict:
    """Strict seven-member/source/type/actual-PCM/full-replay check; never repair."""
    _validate_fixed_contract()
    directory = Path(directory)
    expected_members = set(FILE_NAMES.values()) | {"MANIFEST.json"}
    validate_asset_directory(directory, expected_members, check=True)

    manifest = strict_json_loads((directory / "MANIFEST.json").read_bytes())
    if not isinstance(manifest, dict) or manifest.get("source_sha256") != source_digests():
        raise ValueError("wearable audio source set or SHA differs")

    files = manifest.get("files")
    samples = manifest.get("samples")
    if (
        not isinstance(files, dict)
        or set(files) != set(FILE_NAMES.values())
        or not isinstance(samples, dict)
        or set(samples) != set(FILE_NAMES)
        or any(not isinstance(item, dict) for item in (*files.values(), *samples.values()))
    ):
        raise ValueError("wearable audio manifest member records differ")

    decoded: dict[str, np.ndarray] = {}
    blobs: dict[str, bytes] = {}
    for key, filename in FILE_NAMES.items():
        blob = (directory / filename).read_bytes()
        blobs[filename] = blob
        channels = _FIXED_CHANNELS[key]
        try:
            with wave.open(io.BytesIO(blob)) as reader:
                actual = (
                    reader.getframerate(),
                    reader.getnchannels(),
                    reader.getnframes(),
                    reader.getsampwidth(),
                    reader.getcomptype(),
                )
                if actual != (SAMPLE_RATE, channels, SAMPLES, 2, "NONE"):
                    raise ValueError("wearable audio PCM format differs: " + filename)
            _, decoded[key] = read_pcm16(blob)
        except (wave.Error, EOFError) as exc:
            raise ValueError("invalid wearable audio WAV: " + filename) from exc
        record = {
            "channels": channels,
            "samples_per_channel": SAMPLES,
            "sample_rate_hz": SAMPLE_RATE,
            "sha256": _sha(blob),
        }
        if not same_metadata(files[filename], record) or samples[key].get("file") != filename:
            raise ValueError("wearable audio WAV metadata/SHA differs: " + filename)

    pcm_ref = decoded["wearable_clean_ref"]
    for key, filename in FILE_NAMES.items():
        actual_pcm_score = measure_signal(decoded[key], key, reference=pcm_ref, pcm=True)
        if not same_metadata(samples[key].get("pcm_measurements"), actual_pcm_score):
            raise ValueError("wearable audio PCM scores differ: " + filename)

    expected, contents = prepare_assets()
    if not same_metadata(manifest, expected):
        raise ValueError("wearable audio metadata/model replay differs")
    if any((directory / name).read_bytes() != blob for name, blob in contents.items()):
        raise ValueError("wearable audio byte replay differs from current source/environment")
    return manifest


def generate_assets(directory: Path | str = DEFAULT_OUTPUT, *, check: bool = False) -> dict:
    if type(check) is not bool:
        raise ValueError("check must be bool")
    _validate_fixed_contract()
    directory = Path(directory)
    if check:
        return check_assets(directory)
    expected_members = set(FILE_NAMES.values()) | {"MANIFEST.json"}
    validate_asset_directory(directory, expected_members, check=False)
    manifest, contents = prepare_assets()
    directory.mkdir(parents=True, exist_ok=True)
    validate_asset_directory(directory, expected_members, check=False)
    for filename, blob in contents.items():
        (directory / filename).write_bytes(blob)
    return check_assets(directory)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    result = generate_assets(args.output_dir, check=args.check)
    print(
        json.dumps(
            {
                "status": "checked" if args.check else "generated",
                "files": list(result["files"]),
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
