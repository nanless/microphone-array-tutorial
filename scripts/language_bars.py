#!/usr/bin/env python3
"""Maintain the GitHub-visible 中文 / English switch line at the top of every bilingual Markdown pair.

GitHub cannot run JavaScript inside Markdown, so each document carries one plain line that starts with
``🌐 `` and links to its counterpart in the other language.  The line is *repository navigation only*:
``build_markdown_helpers.shield_math`` removes it before the website and the PDF are rendered (they have
their own language pills), so no published page or PDF text changes.

Usage:
    python scripts/language_bars.py           # insert / refresh every bar
    python scripts/language_bars.py --check   # read-only verification, non-zero exit on drift
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parent.parent
MARK = "🌐"
BAR_LINE = re.compile(r"^🌐 [^\n]*\n?", re.M)

CH = ROOT / "chapters"
CH_EN = ROOT / "chapters_en"
RS = ROOT / "codes/chapters/ch00/research"
RS_EN = ROOT / "codes/chapters/ch00/research_en"

_SUFFIXED = (
    "AGENTS.md",
    "README.md",
    "scripts/README.md",
    "codes/chapters/README.md",
    "codes/chapters/ch00/COVERAGE.md",
    "codes/chapters/ch00/THIRD_PARTY.md",
    "codes/chapters/ch00/upstream/README.md",
    "scripts/vendor/mathjax-3.2.2/README.md",
    "assets/readme/lab/README.md",
)


def pairs() -> list[tuple[Path, Path]]:
    """Return (Chinese, English) document pairs in a stable order."""
    result = [(ROOT / zh, ROOT / zh.replace(".md", "_EN.md")) for zh in _SUFFIXED]
    result += [(path, CH_EN / path.name) for path in sorted(CH.glob("*.md"))]
    result += [(path, RS_EN / path.name) for path in sorted(RS.glob("*.md"))]
    return result


def _relative(target: Path, start: Path) -> str:
    return os.path.relpath(target, start.parent).replace(os.sep, "/")


def bar_for(path: Path, other: Path, *, english: bool) -> str:
    link = _relative(other, path)
    if english:
        return f"{MARK} [简体中文]({link}) ｜ **English**"
    return f"{MARK} **简体中文** ｜ [English]({link})"


def strip_language_bar(text: str) -> str:
    """Remove the repository-only switch line (first occurrence near the top)."""
    match = BAR_LINE.search(text)
    if not match or match.start() > 600:
        return text
    end = match.end()
    blank = re.match(r"[ \t]*\n", text[end:])
    if blank:
        end += blank.end()
    return text[:match.start()] + text[end:]


def with_bar(text: str, bar: str) -> str:
    """Return text carrying exactly one up-to-date bar (after a leading H1, otherwise first line)."""
    text = strip_language_bar(text)
    first, sep, rest = text.partition("\n")
    if first.startswith("# "):
        return f"{first}\n\n{bar}\n\n{rest.lstrip(chr(10))}" if sep else f"{first}\n\n{bar}\n"
    return f"{bar}\n\n{text}"


def sync(check: bool) -> list[str]:
    problems = []
    for zh, en in pairs():
        for path, other, english in ((zh, en, False), (en, zh, True)):
            if not path.is_file() or not other.is_file():
                problems.append(f"missing pair member: {path.relative_to(ROOT)} / {other.relative_to(ROOT)}")
                continue
            old = path.read_text(encoding="utf-8")
            new = with_bar(old, bar_for(path, other, english=english))
            if new != old:
                if check:
                    problems.append(f"language bar missing or stale: {path.relative_to(ROOT)}")
                else:
                    path.write_text(new, encoding="utf-8")
    return problems


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--check", action="store_true", help="verify only; never write")
    args = parser.parse_args(argv)
    problems = sync(args.check)
    for item in problems:
        print(item, file=sys.stderr)
    if problems:
        return 1
    count = sum(2 for _ in pairs())
    print(f"language bars {'verified' if args.check else 'synchronized'}: {count} documents")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
