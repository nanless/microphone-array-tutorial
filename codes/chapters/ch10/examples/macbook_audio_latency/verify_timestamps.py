"""Recheck the archived offsets and callback timestamps without recording audio."""

import csv
import json
from collections import Counter
from decimal import Decimal
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parent
RATE = 48_000
CONFIG = {"null": (128, 8, 1988), "base128": (128, 20, 4688),
          "base480": (480, 20, 1250), "paired": (128, 32, 7388)}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load_trace(name, frames, count):
    with (ROOT / f"{name}.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    require(len(rows) == count, f"{name}: callback count mismatch")
    for index, row in enumerate(rows):
        require((int(row["start"]), int(row["frames"])) == (index * frames, frames),
                f"{name}: callback range mismatch")
    return rows


def sample_time(rows, sample, field, frames):
    callback, offset = divmod(sample, frames)
    return Decimal(rows[callback][field]) + Decimal(offset) / RATE


def main():
    report = json.loads((ROOT / "results.json").read_text(encoding="utf-8"), parse_float=Decimal)
    traces = {}
    for name, (frames, trials, callbacks) in CONFIG.items():
        traces[name] = load_trace(name, frames, callbacks)
        require(len(report[name]["trials"]) == trials, f"{name}: trial count mismatch")
        flags = Counter(row["flags"] for row in traces[name])
        require(flags == report[name]["callback_flags_histogram"], f"{name}: callback flags differ")
        for row in report[name]["trials"]:
            lag = Decimal(row["first_offset_samples"]) * 1000 / RATE
            require(abs(lag - row["first_offset_ms"]) < Decimal("1e-10"),
                    f"{name}: sample-to-time mismatch")

    by_delay = {60: [], 100: []}
    paired_offsets = []
    for index, row in enumerate(report["paired"]["trials"]):
        delay_ms = (60, 100)[index % 2]
        require((row["trial"], row["programmed_delay_ms"]) == (index + 1, delay_ms),
                "unexpected pairing order")
        incoming = 24_000 + index * 28_800 + 4_100
        outgoing = incoming + delay_ms * RATE // 1000
        adc = sample_time(traces["paired"], incoming, "adc_s", 128)
        dac = sample_time(traces["paired"], outgoing, "dac_s", 128)
        delay = (dac - adc) * 1000
        # The recorded analyzer used float64 on timestamps over 10^6 s.
        require(abs(delay - row["host_tagged_adc_to_dac_ms"]) < Decimal("0.000001"),
                f"paired trial {index + 1}: timestamp mismatch")
        by_delay[delay_ms].append(delay)
        echo = delay_ms + row["echo_return_offset_ms"]
        require(abs(echo - row["first_to_forwarded_echo_ms"]) < Decimal("1e-10"),
                f"paired trial {index + 1}: stored offset mismatch")
        paired_offsets.append(echo)

    for name in ("base128", "base480"):
        offset = median(row["first_offset_ms"] for row in report[name]["trials"])
        print(f"{name}: stored acoustic offset {offset:.5f} ms")
    for delay_ms, values in by_delay.items():
        print(f"{delay_ms} ms buffer: recomputed timestamp median {median(values):.6f} ms")
    differences = [b - a for a, b in zip(paired_offsets[::2], paired_offsets[1::2])]
    print(f"Stored acoustic differences: {min(differences):.5f}–{max(differences):.5f} ms")
    print("32 timestamp pairs verified. Correlation peaks are archived values.")


if __name__ == "__main__":
    main()
