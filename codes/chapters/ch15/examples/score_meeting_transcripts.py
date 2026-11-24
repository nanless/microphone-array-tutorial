"""Compare pre-normalized SegLST transcriptions with the same session references.

Each segment has session_id, speaker, start_time and end_time in seconds, and
words as a whitespace-delimited string. Install MeetEval 0.4.3 and simplejson
in a separate environment. The script reads local files and prints scores.
"""

from __future__ import annotations

import argparse
from decimal import Decimal
import json
import math
from pathlib import Path
from statistics import fmean


def _summarize(by_session: dict) -> dict:
    sessions = list(by_session.values())
    reference_words = sum(result.length for result in sessions)
    errors = sum(result.errors for result in sessions)
    if any(result.length == 0 for result in sessions):
        raise ValueError("Every scored session must contain reference words")
    return {
        "session_mean": fmean(result.error_rate for result in sessions),
        "word_weighted": errors / reference_words,
        "reference_words": int(reference_words),
        "errors": int(errors),
        "insertions": int(sum(result.insertions for result in sessions)),
        "deletions": int(sum(result.deletions for result in sessions)),
        "substitutions": int(sum(result.substitutions for result in sessions)),
    }


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, nargs="+", required=True)
    parser.add_argument(
        "--hypothesis",
        nargs="+",
        action="append",
        required=True,
        metavar="VALUE",
        help="System name followed by one or more SegLST paths; repeat per system",
    )
    parser.add_argument("--collar", type=float, default=5.0, help="Seconds")
    args = parser.parse_args(argv)
    if not math.isfinite(args.collar) or args.collar < 0:
        parser.error("--collar must be finite and nonnegative")

    try:
        import simplejson  # Required by MeetEval's SegLST reader.
        from meeteval.io import load
        from meeteval.wer import cpwer, tcpwer
    except ImportError:
        parser.error("Install meeteval==0.4.3 and simplejson==3.20.2 in the scoring environment")

    reference = load(args.reference, format="seglst")
    reference_sessions = set(reference.groupby("session_id"))
    if not reference_sessions:
        parser.error("The reference has no sessions")

    scores = {}
    for name, *paths in args.hypothesis:
        if not paths or name in scores:
            parser.error("Each hypothesis needs a distinct name and at least one file")
        hypothesis = load([Path(path) for path in paths], format="seglst")
        hypothesis_sessions = set(hypothesis.groupby("session_id"))
        if reference_sessions != hypothesis_sessions:
            missing = sorted(reference_sessions - hypothesis_sessions)
            extra = sorted(hypothesis_sessions - reference_sessions)
            parser.error(f"{name} session mismatch: missing={missing}, extra={extra}")
        try:
            scores[name] = {
                "cpwer": _summarize(cpwer(reference, hypothesis)),
                "tcpwer": _summarize(
                    tcpwer(reference, hypothesis, collar=Decimal(str(args.collar)))
                ),
            }
        except ValueError as exc:
            parser.error(f"{name}: {exc}")

    print(json.dumps({
        "sessions": len(reference_sessions),
        "collar_seconds": args.collar,
        "scores": scores,
    }, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
