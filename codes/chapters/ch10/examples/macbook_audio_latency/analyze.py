"""Analyze the PortAudio trials; writes metrics only, no microphone waveforms."""
import json
import sys
from pathlib import Path

import numpy as np

RATE = 48_000
FIRST = 24_000
PERIOD = 28_800
GATE_START = 4_100
GATE_LENGTH = 2_700
CONDITIONS = {
    "null": (8, 0, 128),
    "base128": (20, 0, 128),
    "paired": (32, 2_880, 128),
    "base480": (20, 0, 480),
}


def match(reference, samples, index, lo, hi):
    reference = reference - np.mean(reference)
    n = len(reference)
    window = samples[index + lo:index + hi + n]
    assert len(window) == hi - lo + n
    squares = np.r_[0., np.cumsum(window**2)]
    sums = np.r_[0., np.cumsum(window)]
    energy = squares[n:] - squares[:-n] - (sums[n:] - sums[:-n])**2 / n
    power = (reference @ reference) * energy
    correlation = np.correlate(window, reference, "valid") / np.sqrt(np.maximum(power, 1e-15))
    peak = int(np.argmax(np.abs(correlation)))
    return lo + peak, float(correlation[peak])


def describe(values):
    return {
        "count": len(values),
        "median": float(np.median(values)),
        "min": float(min(values)),
        "max": float(max(values)),
    }


def inspect(directory, name, probes, delay, requested_frames):
    prefix = directory / name
    rx = np.fromfile(prefix.with_suffix(".input.f32"), dtype="<f4").astype(float)
    tx = np.fromfile(prefix.with_suffix(".output.f32"), dtype="<f4").astype(float)
    probe = np.fromfile(prefix.with_suffix(".probe.f32"), dtype="<f4").astype(float)
    t = np.genfromtxt(prefix.with_suffix(".csv"), delimiter=",", names=True)
    expected_samples = FIRST + PERIOD * probes
    assert len(rx) == len(tx) == expected_samples
    assert len(probe) == 960
    assert all(np.isfinite(data).all() for data in (rx, tx, probe))
    assert t["start"][0] == 0 and np.all(np.diff(t["start"]) == t["frames"][:-1])
    assert (t["start"][-1] + t["frames"][-1]) >= len(rx)
    assert np.all(np.diff(t["adc_s"]) > 0) and np.all(np.diff(t["dac_s"]) > 0)
    device = dict(line.split("=", 1) for line in prefix.with_suffix(".device.txt").read_text().splitlines())
    if name == "null":
        assert not np.any(tx)

    def timestamp(field, index):
        cb = np.searchsorted(t["start"], index, side="right") - 1
        return t[field][cb] + (index - t["start"][cb]) / RATE

    rows = []
    for trial in range(probes):
        per_trial_delay = delay + (1_920 if delay and trial % 2 else 0)
        start = FIRST + PERIOD * trial
        if name != "null":
            np.testing.assert_array_equal(tx[start:start+len(probe)], probe)
        lag, corr = match(probe, rx, start, 0, 11_040)
        row = {
            "trial": trial + 1,
            "first_correlation": corr,
            "first_offset_samples": lag,
            "first_offset_ms": lag * 1000 / RATE,
        }
        if delay:
            gate = start + GATE_START
            replay = gate + per_trial_delay
            replay_samples = tx[replay:replay + GATE_LENGTH]
            expected = np.clip(rx[gate:gate + GATE_LENGTH].astype(np.float32) * np.float32(0.28),
                               np.float32(-0.14), np.float32(0.14))
            np.testing.assert_array_equal(replay_samples, expected)
            if not GATE_START <= lag <= GATE_START + GATE_LENGTH - len(probe):
                row["whole_first_probe_within_gate"] = False
            else:
                row["whole_first_probe_within_gate"] = True
            lag2, corr2 = match(replay_samples, rx, replay, 3_600, 7_200)
            row.update({
                "echo_correlation": corr2,
                "echo_return_offset_ms": lag2 * 1000 / RATE,
                "first_to_forwarded_echo_ms": (per_trial_delay + lag2) * 1000 / RATE,
                "programmed_delay_ms": per_trial_delay * 1000 / RATE,
                "host_tagged_adc_to_dac_ms": 1000 * (timestamp("dac_s", replay) - timestamp("adc_s", gate)),
                "forwarding_clip_count": int(np.count_nonzero(np.abs(replay_samples) >= np.float32(0.14))),
            })
        rows.append(row)
    series = {"absolute_first_correlation": describe([abs(row["first_correlation"]) for row in rows])}
    if name != "null":
        series["first_offset_ms"] = describe([row["first_offset_ms"] for row in rows])
    if delay:
        for key in ("echo_correlation", "echo_return_offset_ms", "first_to_forwarded_echo_ms", "host_tagged_adc_to_dac_ms"):
            series[key] = describe([row[key] for row in rows])
    flags, flag_counts = np.unique(t["flags"].astype(int), return_counts=True)
    paired = {}
    if delay:
        keys = ("first_offset_ms", "echo_correlation", "echo_return_offset_ms",
                "first_to_forwarded_echo_ms", "host_tagged_adc_to_dac_ms")
        for value in (60, 100):
            subgroup = [row for row in rows if row["programmed_delay_ms"] == value]
            paired[str(value)] = {key: describe([row[key] for row in subgroup]) for key in keys}
        paired["100_minus_60"] = {
            key: describe([b[key] - a[key] for a, b in zip(rows[::2], rows[1::2])])
            for key in keys
        }
    return {
        "alternating_delay_ms": paired,
        "device": device,
        "requested_callback_frames": requested_frames,
        "observed_callback_frames": np.unique(t["frames"].astype(int)).tolist(),
        "callback_flags_histogram": dict(zip([str(flag) for flag in flags], flag_counts.tolist())),
        "same_index_host_adc_to_dac_ms": describe(1000 * (t["dac_s"] - t["adc_s"])),
        "programmed_delay_ms": delay * 1000 / RATE,
        "summary": series,
        "trials": rows,
    }


if __name__ == "__main__":
    directory = Path(sys.argv[1])
    results = {name: inspect(directory, name, *args) for name, args in CONDITIONS.items()}
    output = directory / "results.json"
    output.write_text(json.dumps(results, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    for name, result in results.items():
        print(name, "flags", result["callback_flags_histogram"],
              "same index host median", result["same_index_host_adc_to_dac_ms"]["median"])
        for metric, summary in result["summary"].items():
            print(" ", metric, summary)
        if result["alternating_delay_ms"]:
            print(" ", "paired summary", result["alternating_delay_ms"])
            print(" ", "whole-probe gates", sum(row["whole_first_probe_within_gate"] for row in result["trials"]),
                  "clipped samples", sum(row["forwarding_clip_count"] for row in result["trials"]))
