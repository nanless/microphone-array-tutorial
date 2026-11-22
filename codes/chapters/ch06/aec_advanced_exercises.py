"""Exact, dimensionless answer key for E06-08..E06-10.

The cases are deliberately tiny. They verify algebra and model boundaries,
not speech quality, detector performance, or industrial AEC throughput.
Importing this module performs no experiment or file access.
"""
from __future__ import annotations

if __name__ == "__main__" and not __package__:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import json


def run_exercises() -> dict:
    # E06-08: IPNLMS pure-proportionate zero-start deadlock at kappa=1.
    # E06-09: exponential forgetting measures UPDATE count, not physical seconds.
    old_weight = 0.99 ** 100
    # E06-10: the fourth physical tap may reach the block two steps back.
    physical_taps = 4
    block_taps = physical_taps // 2 + 1

    return {
        "E06-08": {
            "pure_proportionate_zero_start_update": [0.0, 0.0],
            "reason": "all g_l vanish at kappa=1 and zero weights",
        },
        "E06-09": {
            "old_observation_weight_after_100_updates": old_weight,
            "elapsed_seconds_if_16khz_sample_updates": 100 / 16000,
            "elapsed_seconds_if_10ms_block_updates": 100 * .01,
        },
        "E06-10": {
            "physical_taps": physical_taps,
            "minimum_block_taps": block_taps,
            "unconstrained_crossband_weights": 4 * block_taps,
            "minimum_pairing_wait_samples": 1,
            "scope": "two-band Haar block FIR, not an industrial cost benchmark",
        },
    }


if __name__ == "__main__":
    print(json.dumps(run_exercises(), ensure_ascii=False, indent=2, allow_nan=False))
