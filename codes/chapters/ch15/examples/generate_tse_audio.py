"""Generate six Chapter 15 near-field distance-bubble WAVs; --check replays strictly read-only."""

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
from codes.chapters.ch15.core.tse_audio import (
    CHANNELS_BY_KEY,
    FILE_NAMES,
    LIMITS,
    OUTPUT_KEYS,
    SAMPLE_RATE,
    SAMPLES,
    STEADY_SAMPLES,
    STEADY_START,
    STEADY_STOP,
    measure_signal,
    parameters,
    run_experiment,
)

ROOT = Path(__file__).resolve().parents[4]
DEFAULT_OUTPUT = ROOT / "codes/chapters/ch15/tse_audio"
SOURCE_PATHS = (
    "codes/chapters/ch15/core/distance_bubble.py",
    "codes/chapters/ch15/core/tse_audio.py",
    "codes/chapters/ch15/examples/generate_tse_audio.py",
    "codes/chapters/ch02/core/conventions.py",
    "codes/chapters/ch00/core/audio_samples.py",
    "codes/chapters/ch00/io_contracts.py",
)

_FIXED_FILES = {
    "near_target_ref": "tse_near_target_ref.wav",
    "far_interferer_ref": "tse_far_interferer_ref.wav",
    "array_3ch": "tse_array_3ch.wav",
    "farfield_dsb": "tse_farfield_dsb.wav",
    "bubble_d10cm": "tse_bubble_d10cm.wav",
    "bubble_d2cm_illcond": "tse_bubble_d2cm_illcond.wav",
}
MEMBERS = frozenset(_FIXED_FILES.values()) | {"MANIFEST.json"}

_FIXED_LIMITS = (
    "Deterministic interior-frequency sinusoids (3200/2500/1100/1300/1700 Hz) with exact integer "
    "cycles on [1600, 30400) and known broadside spherical propagation; not natural speech, "
    "a measured room, blind source range estimation, or a subjective listening test. "
    "Digital amplitude squared is not Pa^2; common export gain 1.0 with no per-file normalization, "
    "fitted gain, or post-hoc delay alignment."
)


def _fixed_parameters() -> dict:
    """Independent publication contract for Chapter 15 synthetic audio."""
    return {
        "exercise_ids": ["E15-03", "E15-05", "E15-12"],
        "sample_rate_hz": 16000,
        "samples_per_channel": 32000,
        "fade_samples_each_end": 800,
        "steady_scoring_window": [1600, 30400],
        "steady_scoring_samples": 28800,
        "speed_of_sound_m_per_s": 343.0,
        "target_frequency_hz": 3200.0,
        "interferer_frequency_hz": 2500.0,
        "sensor_floor_frequencies_hz": [1100.0, 1300.0, 1700.0],
        "target_range_m": 0.40,
        "interferer_range_m": 4.00,
        "broadside_angle_deg": 0.0,
        "ula_spacing_good_m": 0.10,
        "ula_spacing_tiny_m": 0.02,
        "tiny_aperture_tikhonov_reg": 2e-4,
        "target_peak_amplitude": 0.16,
        "interferer_peak_amplitude": 0.16,
        "sensor_floor_peak_amplitude": 0.012,
        "target_steady_mean_square": 0.0128,
        "interferer_steady_mean_square": 0.0128,
        "sensor_floor_steady_mean_square": 0.000072,
        "common_export_gain": 1.0,
    }


def _validate_fixed_contract() -> None:
    if (
        type(SAMPLE_RATE) is not int
        or SAMPLE_RATE != 16000
        or type(SAMPLES) is not int
        or SAMPLES != 32000
        or not same_metadata(FILE_NAMES, _FIXED_FILES)
        or not same_metadata(LIMITS, _FIXED_LIMITS)
        or not same_metadata(parameters(), _fixed_parameters())
    ):
        raise ValueError("tse_audio fixed parameters, types, or limits differ")


def _independent_closed_form_taps(d_m: float, reg: float) -> np.ndarray:
    """Closed-form 2x2 Cramers-rule Gram inverse without calling distance_bubble.py."""
    k_tar = 2.0 * np.pi * 3200.0 / 343.0
    delta_r = float(np.hypot(0.40, d_m) - 0.40)
    u = np.exp(-1j * k_tar * delta_r)
    a1 = np.array([u, 1.0 + 0.0j, u], dtype=np.complex128)
    af = np.ones(3, dtype=np.complex128)
    g00 = 3.0 + reg
    g01 = 2.0 * np.conj(u) + 1.0
    g10 = np.conj(g01)
    g11 = 3.0 + reg
    det = g00 * g11 - g01 * g10
    w = a1 * (g11 / det) + af * (-g10 / det)
    om = 2.0 * np.pi * 3200.0 / 16000.0
    h1 = w.imag / np.sin(om)
    h0 = w.real - h1 * np.cos(om)
    return np.column_stack((h0, h1)).astype(np.float64)


def _fixed_waveforms() -> dict[str, np.ndarray]:
    """Independently reconstruct all six declared signals from closed-form expressions."""
    idx = np.arange(32000, dtype=np.float64)
    env = np.minimum(np.minimum(idx / 800.0, (31999.0 - idx) / 800.0), 1.0)
    t_sec = idx / 16000.0

    near = (env * 0.16 * np.cos(2.0 * np.pi * 3200.0 * t_sec))[None, :]
    far = (env * 0.16 * np.cos(2.0 * np.pi * 2500.0 * t_sec))[None, :]
    noises = np.vstack(
        [env * 0.012 * np.sin(2.0 * np.pi * f_n * t_sec) for f_n in (1100.0, 1300.0, 1700.0)]
    )

    def raw_array(d_m: float) -> np.ndarray:
        xs = np.array([-d_m, 0.0, d_m], dtype=np.float64)
        tau1 = (np.hypot(0.40, xs) - 0.40) / 343.0
        tau2 = (np.hypot(4.00, xs) - 4.00) / 343.0
        c1 = env[None, :] * 0.16 * np.cos(2.0 * np.pi * 3200.0 * (t_sec[None, :] - tau1[:, None]))
        c2 = env[None, :] * 0.16 * np.cos(2.0 * np.pi * 2500.0 * (t_sec[None, :] - tau2[:, None]))
        return c1 + c2 + noises

    arr10 = raw_array(0.10)
    arr02 = raw_array(0.02)
    dsb = ((arr10[0] + arr10[1] + arr10[2]) / 3.0)[None, :]

    taps10 = _independent_closed_form_taps(0.10, 0.0)
    bub10 = (taps10[:, 0] @ arr10)[None, :]
    bub10[0, 1:] += taps10[:, 1] @ arr10[:, :-1]

    taps02 = _independent_closed_form_taps(0.02, 2e-4)
    bub02 = (taps02[:, 0] @ arr02)[None, :]
    bub02[0, 1:] += taps02[:, 1] @ arr02[:, :-1]

    return {
        "near_target_ref": near,
        "far_interferer_ref": far,
        "array_3ch": arr10,
        "farfield_dsb": dsb,
        "bubble_d10cm": bub10,
        "bubble_d2cm_illcond": bub02,
    }


def _independent_measurements(x: np.ndarray, key: str, ref_near: np.ndarray, *, pcm: bool = False) -> dict:
    """Direct DFT/energy oracle independent of tse_audio.measure_signal."""
    win = x[:, STEADY_START:STEADY_STOP]
    t_win = np.arange(STEADY_START, STEADY_STOP, dtype=np.float64) / 16000.0
    p_tar = (2.0 / STEADY_SAMPLES) * (win @ np.exp(-1j * 2.0 * np.pi * 3200.0 * t_win))
    p_int = (2.0 / STEADY_SAMPLES) * (win @ np.exp(-1j * 2.0 * np.pi * 2500.0 * t_win))
    out: dict = {
        "steady_window": [1600, 30400],
        "scoring_samples_per_channel": 28800,
        "peak_all_samples": float(np.max(np.abs(x))),
        "mean_square_per_channel": [float(np.sum(r * r) / 28800.0) for r in win],
        "target_phasor_real_imag": [[float(z.real), float(z.imag)] for z in p_tar],
        "interferer_phasor_real_imag": [[float(z.real), float(z.imag)] for z in p_int],
    }
    if key in OUTPUT_KEYS:
        ref_w = ref_near[0, STEADY_START:STEADY_STOP]
        err = win[0] - ref_w
        mse = float(np.sum(err * err) / 28800.0)
        ref_ms = float(np.sum(ref_w * ref_w) / 28800.0)
        out["reference_mse"] = mse
        out["reference_nmse"] = float(mse / ref_ms)
    if pcm:
        ints = np.rint(x * 32768.0).astype(np.int64)[:, STEADY_START:STEADY_STOP]
        e_ch = [int(sum(int(v) * int(v) for v in row)) for row in ints]
        d_ch = 28800 * (32768 ** 2)
        pcm_dict: dict = {
            "pcm_decode_divisor": 32768,
            "integer_squared_sum_E_per_channel": e_ch,
            "integer_denominator_D_per_channel": d_ch,
            "integer_squared_sum_E_all_channels": int(sum(e_ch)),
            "integer_denominator_D_all_channels": int(x.shape[0] * d_ch),
            "mean_square_all_channels": float(sum(e_ch) / (x.shape[0] * d_ch)),
        }
        if key in OUTPUT_KEYS:
            r_int = np.rint(ref_near[0, STEADY_START:STEADY_STOP] * 32768.0).astype(np.int64)
            diff_i = ints[0] - r_int
            e_err = int(sum(int(v) * int(v) for v in diff_i))
            p_ref = int(sum(int(v) * int(v) for v in r_int))
            pcm_dict["integer_reference_error_E"] = e_err
            pcm_dict["integer_reference_power_P"] = p_ref
            pcm_dict["integer_reference_nmse"] = float(e_err / p_ref)
        out["pcm_integer_measurements"] = pcm_dict
    return out


def _numeric_metadata(actual, expected) -> bool:
    """Exact types and keys; allow float64 ULP tolerance only on independent floats."""
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
        return np.isfinite(actual) and bool(np.isclose(actual, expected, rtol=2e-12, atol=1e-13))
    return actual == expected


def _independent_analytic_expectations() -> dict:
    """Closed-form spectral-line oracle using independent Cramer's-rule taps.

    Each interior frequency completes integer cycles in the scoring window;
    squared complex amplitudes therefore add with a factor of one half.
    No waveform projection or imported model expectation is used here.
    """
    source_power = float(0.16 ** 2 / 2)
    noise_power = float(0.012 ** 2 / 2)
    expected = {
        "near_target_ref": {"mean_square_per_channel": [source_power]},
        "far_interferer_ref": {"mean_square_per_channel": [source_power]},
        "array_3ch": {"mean_square_per_channel": [2 * source_power + noise_power] * 3},
    }
    controls = (
        ("farfield_dsb", 0.10, np.column_stack((np.ones(3) / 3, np.zeros(3)))),
        ("bubble_d10cm", 0.10, _independent_closed_form_taps(0.10, 0.0)),
        ("bubble_d2cm_illcond", 0.02, _independent_closed_form_taps(0.02, 2e-4)),
    )
    for key, spacing, taps in controls:
        def response(frequency, distance):
            total = 0j
            for position, (h0, h1) in zip((-spacing, 0.0, spacing), taps):
                path = (float(np.hypot(distance, position)) - distance) / 343.0
                total += (h0 + h1 * np.exp(-2j * np.pi * frequency / 16000)) * np.exp(-2j * np.pi * frequency * path)
            return total
        target = response(3200, 0.4)
        interference = response(2500, 4.0)
        sensor = float(sum(noise_power * abs(h0 + h1 * np.exp(-2j * np.pi * f / 16000)) ** 2
                           for f, (h0, h1) in zip((1100, 1300, 1700), taps)))
        target_power = float(source_power * abs(target) ** 2)
        distortion = float(source_power * abs(target - 1) ** 2)
        leakage = float(source_power * abs(interference) ** 2)
        mse = distortion + leakage + sensor
        expected[key] = {
            "mean_square_per_channel": [target_power + leakage + sensor],
            "target_output_power": target_power,
            "target_distortion_power": distortion,
            "interferer_leakage_power": leakage,
            "sensor_noise_output_power": sensor,
            "reference_mse": mse,
            "reference_nmse": float(mse / source_power),
            "interferer_rejection_db": float(-10 * np.log10(max(abs(interference) ** 2, 1e-30))),
        }
    return expected


def _validate_model(report: dict, arrays: dict[str, np.ndarray]) -> None:
    fixed = _fixed_waveforms()
    if (
        type(report) is not dict
        or not same_metadata(report.get("parameters"), _fixed_parameters())
        or not same_metadata(report.get("limits"), _FIXED_LIMITS)
    ):
        raise ValueError("tse_audio model fixed metadata/types/limits differ")
    if not _numeric_metadata(report.get("analytic_expectations"), _independent_analytic_expectations()):
        raise ValueError("tse_audio independent analytic expectations differ")
    for row in report["analytic_expectations"].values():
        for field, value in row.items():
            if any(token in field for token in ("power", "mse", "mean_square")):
                values = value if isinstance(value, list) else [value]
                if any(v < 0.0 for v in values):
                    raise ValueError("independent analytic powers/errors must be non-negative")
    if type(arrays) is not dict or arrays.keys() != fixed.keys():
        raise ValueError("tse_audio fixture members differ")
    for key, expected in fixed.items():
        arr = arrays[key]
        if (
            type(arr) is not np.ndarray
            or arr.dtype != np.dtype("float64")
            or arr.shape != expected.shape
            or not np.all(np.isfinite(arr))
            or not np.allclose(arr, expected, rtol=0.0, atol=1e-13)
        ):
            raise ValueError(f"tse_audio fixed waveform/type/shape differs: {key}")
    expected_float = {
        k: _independent_measurements(x, k, fixed["near_target_ref"], pcm=False)
        for k, x in fixed.items()
    }
    if not _numeric_metadata(report.get("float_measurements"), expected_float):
        raise ValueError("tse_audio independent float measurements differ")


def _sha(blob: bytes) -> str:
    return hashlib.sha256(blob).hexdigest()


def source_digests() -> dict[str, str]:
    digests: dict[str, str] = {}
    for path in SOURCE_PATHS:
        source = validate_parent_chain(ROOT / path)
        if not source.is_file():
            raise ValueError(f"missing ordinary audio source: {path}")
        digests[path] = _sha(source.read_bytes())
    return digests


def prepare_assets() -> tuple[dict, dict[str, bytes]]:
    """Build all six WAV byte blobs and MANIFEST.json in memory after independent checks."""
    _validate_fixed_contract()
    report, arrays = run_experiment()
    _validate_model(report, arrays)

    contents: dict[str, bytes] = {}
    files: dict[str, dict] = {}
    decoded: dict[str, np.ndarray] = {}
    for key, arr in arrays.items():
        blob = pcm16_bytes(arr, SAMPLE_RATE)
        rate, pcm = read_pcm16(blob)
        channels = CHANNELS_BY_KEY[key]
        with wave.open(io.BytesIO(blob), "rb") as reader:
            if (
                reader.getframerate(),
                reader.getnchannels(),
                reader.getnframes(),
                reader.getsampwidth(),
                reader.getcomptype(),
            ) != (16000, channels, 32000, 2, "NONE"):
                raise ValueError(f"tse_audio in-memory WAV format differs: {key}")
            integers = np.frombuffer(reader.readframes(32000), dtype="<i2").reshape(-1, channels).T
        if (
            type(rate) is not int
            or rate != 16000
            or not np.array_equal(integers, np.rint(arr * 32768.0))
            or type(pcm) is not np.ndarray
            or pcm.dtype != np.dtype("float64")
            or not np.array_equal(pcm, integers.astype(np.float64) / 32768.0)
        ):
            raise ValueError(f"tse_audio independent PCM encoding/decoding differs: {key}")
        filename = FILE_NAMES[key]
        contents[filename] = blob
        decoded[key] = pcm
        files[filename] = {
            "channels": int(channels),
            "samples_per_channel": SAMPLES,
            "sample_rate_hz": rate,
            "sha256": _sha(blob),
        }

    samples: dict[str, dict] = {}
    for key, arr in arrays.items():
        filename = FILE_NAMES[key]
        pcm = decoded[key]
        pcm_scores = measure_signal(pcm, key, reference=decoded["near_target_ref"][0], pcm=True)
        indep_pcm = _independent_measurements(pcm, key, decoded["near_target_ref"], pcm=True)
        if not _numeric_metadata(pcm_scores, indep_pcm):
            raise ValueError(f"tse_audio independent PCM metadata/score differs: {key}")
        samples[key] = {
            "file": filename,
            "analytic": report["analytic_expectations"][key],
            "float_measurements": report["float_measurements"][key],
            "pcm_measurements": pcm_scores,
            "quantization_max_abs_error": float(np.max(np.abs(arr - pcm))),
        }

    manifest = {
        "schema_version": 1,
        "origin": "original synthetic broadside near-field distance-bubble controls",
        "sample_rate_hz": SAMPLE_RATE,
        "samples_per_channel": SAMPLES,
        "common_export_gain": 1.0,
        "source_sha256": source_digests(),
        "parameters": parameters(),
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
        "pcm_half_step_error_bound": 0.5 / 32768.0,
        "pcm": "signed little-endian PCM16, nearest-even rounding; no dither, saturation, or per-file gain",
        "limits": LIMITS,
        "listening": "Start at low device volume; no autoplay or subjective listening claim.",
    }
    contents["MANIFEST.json"] = (
        json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    ).encode("utf-8")
    return manifest, contents


def check_assets(directory: Path | str = DEFAULT_OUTPUT) -> dict:
    """Strict 7-member, SHA-256, WAV-header, PCM-score, and full-replay read-only check."""
    _validate_fixed_contract()
    directory = Path(directory)
    validate_asset_directory(directory, MEMBERS, check=True)
    manifest = strict_json_loads((directory / "MANIFEST.json").read_bytes())
    if not isinstance(manifest, dict) or manifest.get("source_sha256") != source_digests():
        raise ValueError("tse_audio source set or SHA-256 differs")

    files = manifest.get("files")
    samples = manifest.get("samples")
    if (
        not isinstance(files, dict)
        or set(files) != set(FILE_NAMES.values())
        or not isinstance(samples, dict)
        or set(samples) != set(FILE_NAMES)
        or any(not isinstance(item, dict) for item in (*files.values(), *samples.values()))
    ):
        raise ValueError("tse_audio manifest member records differ")

    decoded: dict[str, np.ndarray] = {}
    raw_blobs: dict[str, bytes] = {}
    for key, filename in FILE_NAMES.items():
        blob = (directory / filename).read_bytes()
        raw_blobs[filename] = blob
        channels = CHANNELS_BY_KEY[key]
        try:
            with wave.open(io.BytesIO(blob), "rb") as reader:
                actual_fmt = (
                    reader.getframerate(),
                    reader.getnchannels(),
                    reader.getnframes(),
                    reader.getsampwidth(),
                    reader.getcomptype(),
                )
                if actual_fmt != (SAMPLE_RATE, channels, SAMPLES, 2, "NONE"):
                    raise ValueError(f"tse_audio PCM format differs: {filename}")
            _, decoded[key] = read_pcm16(blob)
        except (wave.Error, EOFError) as exc:
            raise ValueError(f"invalid tse_audio WAV: {filename}") from exc
        record = {
            "channels": int(channels),
            "samples_per_channel": SAMPLES,
            "sample_rate_hz": SAMPLE_RATE,
            "sha256": _sha(blob),
        }
        if not same_metadata(files[filename], record) or samples[key].get("file") != filename:
            raise ValueError(f"tse_audio WAV metadata/SHA differs: {filename}")

    for key, filename in FILE_NAMES.items():
        actual_pcm_score = measure_signal(
            decoded[key], key, reference=decoded["near_target_ref"][0], pcm=True
        )
        if not same_metadata(samples[key].get("pcm_measurements"), actual_pcm_score):
            raise ValueError(f"tse_audio PCM scores differ: {filename}")

    expected, contents = prepare_assets()
    if not same_metadata(manifest, expected):
        raise ValueError("tse_audio metadata/model replay differs")
    if any((directory / name).read_bytes() != blob for name, blob in contents.items()):
        raise ValueError("tse_audio byte replay differs from current source/environment")
    return manifest


def generate_assets(directory: Path | str = DEFAULT_OUTPUT, *, check: bool = False) -> dict:
    """Generate or strictly check the 7 ordinary files in ``tse_audio/``."""
    if type(check) is not bool:
        raise ValueError("check must be bool")
    _validate_fixed_contract()
    directory = Path(directory)
    if check:
        return check_assets(directory)
    validate_asset_directory(directory, MEMBERS, check=False)
    _manifest, contents = prepare_assets()
    directory.mkdir(parents=True, exist_ok=True)
    validate_asset_directory(directory, MEMBERS, check=False)
    for filename, blob in contents.items():
        (directory / filename).write_bytes(blob)
    return check_assets(directory)


def main(argv: list[str] | None = None) -> int:
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
