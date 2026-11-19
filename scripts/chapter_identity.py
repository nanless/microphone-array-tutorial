"""Published topic identities survive the 12--17 chapter-number permutations.

Display numbers and code directories follow the current book. Hidden PDF names
retain their original topic, and the eight legacy HTML routes remain aliases.
The fixed snapshot records real previous headings, not ordinal guesses.
"""
from __future__ import annotations

import json
from pathlib import Path
import re

def _build_topic_tables(numbering: dict) -> tuple[dict[str, dict], dict[str, str]]:
    topics: dict[str, dict] = {}
    pdf_ids: set[int] = set()
    for row in numbering["chapters"]:
        source = row["source"]
        if source in topics:
            raise ValueError(f"Duplicate chapter source in numbering table: {source}")
        pid = int(row["pdf_identity"])
        if pid in pdf_ids:
            raise ValueError(f"Duplicate pdf_identity in numbering table: {pid}")
        pdf_ids.add(pid)
        topics[source] = row

    raw_routes: dict[str, str] = {}
    for row in topics.values():
        target = row["source"]
        for old in row.get("legacy_sources", [row["original_source"]]):
            if old == target:
                continue
            if old in raw_routes and raw_routes[old] != target:
                raise ValueError(f"Conflicting legacy route {old}: {raw_routes[old]} vs {target}")
            raw_routes[old] = target

    resolved: dict[str, str] = {}
    for old in raw_routes:
        seen = {old}
        cur = raw_routes[old]
        while cur in raw_routes:
            if cur in seen:
                raise ValueError(f"Cyclic legacy chapter route involving {old}")
            seen.add(cur)
            cur = raw_routes[cur]
        resolved[old] = cur
    return topics, resolved


NUMBERING_PATH = Path(__file__).with_name("chapter_numbering.json")
NUMBERING = json.loads(NUMBERING_PATH.read_text(encoding="utf-8"))
TOPICS, ORIGINAL_TO_CURRENT = _build_topic_tables(NUMBERING)
LEGACY_HTML_ROUTES = {
    old.replace(".md", ".html"): new.replace(".md", ".html")
    for old, new in ORIGINAL_TO_CURRENT.items()
}


def canonical_source(name: str) -> str:
    return ORIGINAL_TO_CURRENT.get(name, name)


def historical_source(name: str) -> str:
    row = TOPICS.get(canonical_source(name))
    return row["original_source"] if row else name


def pdf_topic_number(name: str) -> int:
    row = TOPICS.get(canonical_source(name))
    if row:
        return row["pdf_identity"]
    match = re.fullmatch(r"(\d{2})_[^/]+\.md", name)
    if match is None:
        raise ValueError("chapter source requires a two-digit identity")
    return int(match.group(1))


def numbering_heading_aliases(name: str) -> dict[str, str]:
    return NUMBERING["heading_aliases"].get(canonical_source(name), {})


def exercise_alias_records(name: str) -> dict[str, dict[str, object]]:
    return NUMBERING["exercise_aliases"].get(canonical_source(name), {})


def canonical_exercise(name: str, fragment: str) -> str:
    records = exercise_alias_records(name)
    matches: list[str] = []
    for current, row in records.items():
        if fragment in row.get("legacy_ids", [row["original_id"]]) and fragment != current:
            matches.append(current)
    if len(matches) > 1:
        raise ValueError(f"Conflicting exercise alias {fragment} in {name}: {matches}")
    if len(matches) == 1:
        seen = {fragment}
        cur = matches[0]
        while True:
            next_matches = [
                cand for cand, row in records.items()
                if cur in row.get("legacy_ids", [row["original_id"]]) and cand != cur
            ]
            if not next_matches:
                return cur
            if len(next_matches) > 1 or next_matches[0] in seen:
                raise ValueError(f"Cyclic or ambiguous exercise alias chain for {fragment} in {name}")
            seen.add(cur)
            cur = next_matches[0]
    return fragment


def inject_exercise_aliases(html: str, name: str, prefix: str = "") -> str:
    """Attach old topic-local IDs without adding exercises to Markdown counts."""
    records = exercise_alias_records(name)

    def alias(match):
        current = match.group(2)
        row = records.get(current)
        if not row:
            return match.group(0)
        spans = "".join(
            f'<span id="{prefix + old_id}" class="anchor-alias exercise-alias" '
            f'aria-hidden="true"></span>'
            for old_id in row.get("legacy_ids", [row["original_id"]])
            if old_id != current
        )
        return spans + match.group(0)

    return re.sub(
        r'(<a\s+id=["\']' + re.escape(prefix) + r')(e\d{2}-\d{2})(["\']\s*>)',
        alias,
        html,
    )
