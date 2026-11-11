"""E10-16/17 and E11-09: event boundaries, timestamp scope, paired evidence.

Deterministic teaching inputs; no device, model, network or file writes.
Run with ``python -m codes.chapters.ch00.cross_chapter.engineering_boundary_exercises``.
"""

from __future__ import annotations

# Allow the documented direct-file command as well as python -m.
if __name__ == "__main__" and not __package__:
    import sys as _chapter_entry_sys
    from pathlib import Path as _ChapterEntryPath
    _chapter_entry_sys.path.insert(0, str(_ChapterEntryPath(__file__).resolve().parents[4]))

import json
import math
from numbers import Integral, Rational, Real

import numpy as np

from codes.chapters.ch10.core.engineering import simulate_deadline_queue


def callback_timing_ns(input_adc_ns, callback_ns, output_dac_ns):
    """Differences of three ordered timestamps in one known clock domain.

    Inputs are constructed integer nanoseconds, not a claim about device
    accuracy. The first-sample latency requires matching input/output samples;
    the helper cannot establish that correspondence or clock calibration.
    """
    values = (input_adc_ns, callback_ns, output_dac_ns)
    if any(isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral)
           for value in values):
        raise ValueError("timestamps must be integer nanoseconds")
    adc, callback, dac = map(int, values)
    if not adc <= callback <= dac:
        raise ValueError("this example requires adc <= callback <= dac")
    try:
        result = {
            "input_age_ms": (callback - adc) / 1_000_000,
            "output_lead_ms": (dac - callback) / 1_000_000,
            "matched_first_sample_latency_ms": (dac - adc) / 1_000_000,
        }
    except OverflowError as error:
        raise ValueError("timestamp differences exceed finite float64 range") from error
    if not all(math.isfinite(value) for value in result.values()):
        raise ValueError("timestamp differences exceed finite float64 range")
    return result


def variable_callback_timing_demo():
    """Replay three planned markers on a fixed, synthetic 16 kHz duplex trace."""
    sample_rate_hz = 16_000
    sample_period_ns = 62_500
    ns_per_ms = 1_000_000
    # After the overrun, the continuous input index is unknown. Output indices
    # count the application's submitted slots, including untagged filler.
    callbacks = [
        (0, 0, 160, 100, 116, 130, False),
        (160, 160, 80, 110, 124, 140, False),
        (240, 240, 160, 115, 132, 145, False),
        (400, 400, 160, 125, 144, 155, False),
        (None, 560, 80, 137, 154, 165, True),
    ]
    trace = [
        dict(callback_id=callback_id, input_start=input_start, output_start=output_start,
             frame_count=frames, input_adc_ns=adc*ns_per_ms, callback_ns=begin*ns_per_ms,
             output_dac_ns=dac*ns_per_ms, input_overflow=overflow)
        for callback_id, (input_start, output_start, frames, adc, begin, dac, overflow)
        in enumerate(callbacks)
    ]
    # Actual output tags are given by the trace, not inferred from timing.
    delivered = {"A": (1, 20), "B": (3, 0)}
    markers = []
    for label, input_index, ready_ms in (("A", 120, 122), ("B", 200, 133), ("C", 480, 150)):
        capture = next(row for row in trace if row["input_start"] is not None
                       and row["input_start"] <= input_index < row["input_start"] + row["frame_count"])
        input_offset = input_index - capture["input_start"]
        input_time = capture["input_adc_ns"] + input_offset*sample_period_ns
        ready_ns = ready_ms*ns_per_ms
        first_candidate = next(row for row in trace if row["callback_ns"] >= ready_ns)
        marker = dict(label=label, input_index=input_index, input_offset=input_offset,
                      capture_callback=capture["callback_id"], adc_sample_ns=input_time,
                      ready_ns=ready_ns, first_callback_at_or_after_ready=first_candidate["callback_id"],
                      output_callback=None, output_offset=None, output_index=None, dac_sample_ns=None,
                      latency_ns=None)
        output_slot = delivered.get(label)
        if output_slot is None:
            if not first_candidate["input_overflow"]:
                raise AssertionError("this fixture expects to discard the pending marker at an input gap")
            marker["status"] = "unmatched_after_gap"
        else:
            output_callback, offset = output_slot
            playback = trace[output_callback]
            if (playback["callback_ns"] < ready_ns or not 0 <= offset < playback["frame_count"]
                    or any(row["input_overflow"]
                           for row in trace[capture["callback_id"] + 1:output_callback + 1])):
                raise AssertionError("the tagged fixture output must be ready and within a continuous epoch")
            output_time = playback["output_dac_ns"] + offset*sample_period_ns
            marker.update(status="matched", output_callback=output_callback, output_offset=offset,
                          output_index=playback["output_start"] + offset, dac_sample_ns=output_time,
                          latency_ns=output_time - input_time)
        markers.append(marker)

    return {
        "scope": "synthetic callback replay; no audio API or hardware execution",
        "sample_rate_hz": sample_rate_hz,
        "trace": trace,
        "markers": markers,
        "known_gap_sample_count": None,
        "planned_markers": len(markers),
        "matched_markers": sum(marker["status"] == "matched" for marker in markers),
    }


def paired_sign_test_lower_is_better(baseline, candidate):
    """One-sided exact sign test conditional on non-ties, H0: P(win)=1/2.

    Independent paired observations, a prespecified direction, and the stated
    null are prerequisites. Magnitudes and pooled WER are not tested. All ties
    return no p-value; ties remain in the reported total. Equality is exact on
    the supplied observations, not determined by an implicit tolerance.
    """
    def observations(values):
        # Object conversion preserves each Python integer even in mixed lists.
        array = np.asarray(values, dtype=object)
        if array.ndim != 1:
            raise ValueError("observations must be one-dimensional")
        result = []
        for value in array:
            if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
                raise ValueError("observations must be finite real numbers, not booleans")
            if isinstance(value, Integral):
                result.append(int(value))
            elif isinstance(value, Rational):
                # Preserve exact Fraction order, including mixed int/float pairs.
                result.append(value)
            else:
                if not isinstance(value, (float, np.floating)) or (
                        isinstance(value, np.floating) and value.dtype.itemsize > 8):
                    raise ValueError("only integers, rationals and <=64-bit floating observations are supported")
                value = float(value)
                if not math.isfinite(value):
                    raise ValueError("observations must be finite")
                result.append(value)
        return np.asarray(result, dtype=object)
    a, b = observations(baseline), observations(candidate)
    if a.shape != b.shape or not a.size:
        raise ValueError("paired arrays must have equal non-zero length")
    wins = int(np.count_nonzero(b < a))
    losses = int(np.count_nonzero(b > a))
    informative = wins + losses
    # Compare observations without subtracting: opposite large finite values
    # must not overflow merely to determine a sign.
    tail = (sum(math.comb(informative, k) for k in range(wins, informative + 1))
            if informative else None)
    denominator = (1 << informative) if informative else None
    pvalue = tail / denominator if informative else None
    log_pvalue = (math.log(pvalue) if pvalue else
                  math.log(tail) - informative * math.log(2)) if informative else None
    underflow = informative > 0 and pvalue == 0
    return {
        "pairs": int(a.size), "wins": wins, "losses": losses,
        "ties": int(a.size) - informative, "informative_pairs": informative,
        "one_sided_pvalue": pvalue,
        "log_one_sided_pvalue": log_pvalue,
        "pvalue_underflow": underflow,
        "status": ("positive_pvalue_underflow" if underflow else "valid")
                  if informative else "no_non_tied_pairs",
    }


def run_exercises():
    period = 0.1
    adc, callback, dac = 5_000_000_000, 5_015_000_000, 5_045_000_000
    a = np.array([12, 8, 20, 5, 15, 10])
    b = np.array([10, 9, 14, 5, 12, 9])
    return {
        "E10-16": {
            "period_ms": period, "frames": 30,
            "equal_duration": simulate_deadline_queue(np.full(30, period), period, 1),
            "longer_duration": simulate_deadline_queue(np.full(30, 0.11), period, 1),
        },
        "E10-17": {
            "sample_rate_hz": 16000, "frame_samples": 160,
            "input_adc_ns": adc, "callback_ns": callback, "output_dac_ns": dac,
            **callback_timing_ns(adc, callback, dac),
            "compute_ms": 3.0,
            "added_output_block_ms": 10.0,
            "latency_after_one_extra_output_block_ms": 55.0,
            "variable_callback_replay": variable_callback_timing_demo(),
        },
        "E11-09": {
            "reference_words_per_session": [100] * 6,
            "baseline_errors": a.tolist(), "candidate_errors": b.tolist(),
            "baseline_pooled_wer": int(a.sum()) / 600,
            "candidate_pooled_wer": int(b.sum()) / 600,
            "reduction_percentage_points": int(a.sum() - b.sum()) / 6,
            "relative_error_reduction": int(a.sum() - b.sum()) / int(a.sum()),
            **paired_sign_test_lower_is_better(a, b),
        },
    }


if __name__ == "__main__":
    print(json.dumps(run_exercises(), ensure_ascii=False, indent=2, allow_nan=False))
