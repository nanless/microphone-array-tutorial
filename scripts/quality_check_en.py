#!/usr/bin/env python3
"""Comprehensive quality gate for the English Edition, Dual License, and <50 MiB PDFs."""
from __future__ import annotations

from collections import Counter
from html import unescape
from html.parser import HTMLParser
from pathlib import Path
import re
import sys
from urllib.parse import unquote, urlparse
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pypdf import PdfReader
from scripts import build_pdf, build_site_en, quality_check
from scripts.build_markdown_helpers import code_block_ranges, parsed_markdown_headings, shield_math

CHAPTERS_ZH = ROOT / "chapters"
CHAPTERS_EN = ROOT / "chapters_en"
RESEARCH_ZH = ROOT / "codes" / "chapters" / "ch00" / "research"
RESEARCH_EN = ROOT / "codes" / "chapters" / "ch00" / "research_en"
FIGURES_ZH = ROOT / "figures"
FIGURES_EN = ROOT / "figures_en"
SITE_EN = ROOT / "site_en"
DIST = ROOT / "dist"

CJK_RE = re.compile(r"[\u4e00-\u9fff]")
BACKTICK_MATH_RE = re.compile(r"(?<!`)`\$[^`\n]+\$`(?!`)")
TAG_RE = re.compile(r"\\tag\{([^}]+)\}")
EXPLICIT_ID_RE = re.compile(r'<a\s+id=["\']([^"\']+)["\']\s*>')

EXPECTED_EN_DOCS = (
    "README_EN.md",
    "AGENTS_EN.md",
    "scripts/README_EN.md",
    "codes/chapters/README_EN.md",
    "codes/chapters/ch00/COVERAGE_EN.md",
    "codes/chapters/ch00/THIRD_PARTY_EN.md",
    "codes/chapters/ch00/upstream/README_EN.md",
    "scripts/vendor/mathjax-3.2.2/README_EN.md",
    "assets/readme/lab/README_EN.md",
)

EXPECTED_EN_SVGS = (
    "assets/readme/course-cover-en.svg",
    "assets/readme/course-cover-still-en.svg",
    "assets/readme/delay-and-sum-en.svg",
    "assets/readme/listening-room-en.svg",
    "assets/readme/listening-room-still-en.svg",
)


def strip_fenced(md: str) -> str:
    ranges = code_block_ranges(md)
    if not ranges:
        return md
    out = []
    pos = 0
    for s, e in ranges:
        out.append(md[pos:s])
        pos = e
    out.append(md[pos:])
    return "".join(out)


def check_licenses(errors: list[str]) -> None:
    lic = ROOT / "LICENSE"
    lic_docs = ROOT / "LICENSE-DOCS"
    if not lic.is_file() or "MIT License" not in lic.read_text("utf-8"):
        errors.append("Missing or invalid root LICENSE (expected MIT License for original code)")
    if not lic_docs.is_file() or "Attribution-NonCommercial-ShareAlike 4.0" not in lic_docs.read_text("utf-8"):
        errors.append("Missing or invalid LICENSE-DOCS (expected CC BY-NC-SA 4.0 legal text)")


def check_markdown_links(md_text: str, source_path: Path, errors: list[str]) -> None:
    shielded, _repo = shield_math(strip_fenced(md_text))
    prose_no_inline = re.sub(r"`[^`\n]*`", "", shielded)
    for m in re.finditer(r"!?\[[^\]]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)", prose_no_inline):
        raw_href = m.group(1).strip()
        if raw_href.startswith("<") and raw_href.endswith(">"):
            raw_href = raw_href[1:-1]
        parsed = urlparse(raw_href)
        if parsed.scheme or parsed.netloc or not parsed.path:
            continue
        target = (source_path.parent / unquote(parsed.path)).resolve()
        if not target.exists():
            errors.append(f"Broken local link in {source_path.relative_to(ROOT)}: {raw_href}")


def check_paired_corpus(zh_dir: Path, en_dir: Path, label: str, errors: list[str]) -> None:
    zh_files = sorted(p.name for p in zh_dir.glob("*.md"))
    en_files = sorted(p.name for p in en_dir.glob("*.md")) if en_dir.is_dir() else []
    if en_files != zh_files:
        errors.append(
            f"{label} file set mismatch: missing={sorted(set(zh_files) - set(en_files))}, "
            f"extra={sorted(set(en_files) - set(zh_files))}"
        )
        return
    for name in zh_files:
        zh_path = zh_dir / name
        en_path = en_dir / name
        zh_text = zh_path.read_text("utf-8")
        en_text = en_path.read_text("utf-8")
        rel = en_path.relative_to(ROOT)

        zh_heads = [lvl for lvl, _ in parsed_markdown_headings(zh_text)]
        en_heads = [lvl for lvl, _ in parsed_markdown_headings(en_text)]
        if en_heads != zh_heads:
            errors.append(
                f"{rel}: heading hierarchy mismatch vs Chinese ({len(en_heads)} vs {len(zh_heads)} headings)"
            )

        zh_tags = Counter(TAG_RE.findall(strip_fenced(zh_text)))
        en_tags = Counter(TAG_RE.findall(strip_fenced(en_text)))
        if en_tags != zh_tags:
            errors.append(f"{rel}: equation \\tag mismatch vs Chinese: diff={dict((zh_tags - en_tags) | (en_tags - zh_tags))}")

        zh_ids = set(EXPLICIT_ID_RE.findall(strip_fenced(zh_text)))
        en_ids = set(EXPLICIT_ID_RE.findall(strip_fenced(en_text)))
        if en_ids != zh_ids:
            errors.append(f"{rel}: explicit <a id> mismatch: missing={sorted(zh_ids - en_ids)}, extra={sorted(en_ids - zh_ids)}")

        en_prose = strip_fenced(en_text)
        for lno, line in enumerate(en_prose.splitlines(), 1):
            if re.search(r"\\\(|\\\)|\\\[|\\\]", line):
                errors.append(f"{rel}:{lno}: incompatible math delimiter \\( or \\[")
            if BACKTICK_MATH_RE.search(line):
                errors.append(f"{rel}:{lno}: backtick-wrapped math `$...$`")
            cleaned = re.sub(r"\[(?:简体)?中文(?:版)?\]\([^)]+\)", "", line)
            if CJK_RE.search(cleaned):
                errors.append(f"{rel}:{lno}: unexpected CJK characters: {line.strip()[:80]}")

        check_markdown_links(en_text, en_path, errors)


def check_en_repo_docs(errors: list[str]) -> None:
    for rel in EXPECTED_EN_DOCS:
        p = ROOT / rel
        if not p.is_file() or p.stat().st_size == 0:
            errors.append(f"Missing or empty English guide: {rel}")
            continue
        txt = p.read_text("utf-8")
        prose = strip_fenced(txt)
        for lno, line in enumerate(prose.splitlines(), 1):
            if BACKTICK_MATH_RE.search(line):
                errors.append(f"{rel}:{lno}: backtick-wrapped math `$...$`")
            cleaned = re.sub(r"\[(?:简体)?中文(?:版)?\]\([^)]+\)", "", line)
            cleaned = cleaned.replace("EN/中文", "").replace("EN / 中文", "")
            if CJK_RE.search(cleaned):
                errors.append(f"{rel}:{lno}: unexpected CJK characters: {line.strip()[:80]}")
        check_markdown_links(txt, p, errors)


def check_figures_and_svgs(errors: list[str]) -> None:
    zh_figs = sorted(p.name for p in FIGURES_ZH.glob("fig*.png"))
    en_figs = sorted(p.name for p in FIGURES_EN.glob("fig*.png")) if FIGURES_EN.is_dir() else []
    if len(zh_figs) != 92 or en_figs != zh_figs:
        errors.append(
            f"figures_en mismatch: expected 92 matching files, got {len(en_figs)} "
            f"(missing={sorted(set(zh_figs) - set(en_figs))[:8]})"
        )
    else:
        for name in en_figs:
            p = FIGURES_EN / name
            data = p.read_bytes()
            if len(data) < 1024 or not data.startswith(b"\x89PNG\r\n\x1a\n"):
                errors.append(f"Invalid PNG in figures_en/{name}")

    for rel in EXPECTED_EN_SVGS:
        p = ROOT / rel
        if not p.is_file() or p.stat().st_size == 0:
            errors.append(f"Missing English SVG: {rel}")
            continue
        raw = p.read_text("utf-8")
        try:
            ET.fromstring(raw)
        except ET.ParseError as exc:
            errors.append(f"Invalid XML in {rel}: {exc}")
        if CJK_RE.search(raw):
            errors.append(f"Unexpected CJK characters in {rel}")


def check_site_en(errors: list[str]) -> None:
    expected = build_site_en.expected_html_set_en()
    if not SITE_EN.is_dir():
        errors.append("Missing site_en/ directory")
        return
    actual = {p.relative_to(SITE_EN).as_posix() for p in SITE_EN.rglob("*.html")}
    if actual != expected:
        errors.append(f"site_en HTML set mismatch: missing={sorted(expected - actual)}, extra={sorted(actual - expected)}")
        return

    digest = build_site_en.source_digest_en()
    parsed_pages = {}
    for rel in sorted(actual):
        p = SITE_EN / rel
        parser = quality_check.PageParser()
        parser.feed(p.read_text("utf-8"))
        parsed_pages[p.resolve()] = parser
        if parser.html_lang != "en":
            errors.append(f"site_en/{rel}: expected lang='en', got {parser.html_lang!r}")
        if parser.h1_count != 1:
            errors.append(f"site_en/{rel}: expected exactly 1 <h1>, got {parser.h1_count}")
        if parser.source_digest != digest:
            errors.append(f"site_en/{rel}: stale source-digest ({parser.source_digest} != {digest})")
        dups = [k for k, c in Counter(parser.ids).items() if c > 1]
        if dups:
            errors.append(f"site_en/{rel}: duplicate HTML ids: {dups[:5]}")

    for path, parser in parsed_pages.items():
        rel = path.relative_to(SITE_EN.resolve()).as_posix()
        for audio_src in parser.audios:
            aurl = urlparse(audio_src)
            if not aurl.scheme and not aurl.netloc:
                atarget = (path.parent / unquote(aurl.path)).resolve()
                if not atarget.is_file():
                    errors.append(f"site_en/{rel}: missing audio {audio_src}")
        for src, alt in parser.images:
            iurl = urlparse(src)
            if not iurl.scheme and not iurl.netloc:
                itarget = (path.parent / unquote(iurl.path)).resolve()
                if not itarget.is_file():
                    errors.append(f"site_en/{rel}: missing image {src}")
            if not alt or not alt.strip():
                errors.append(f"site_en/{rel}: missing image alt on {src}")
        for href in parser.links:
            url = urlparse(href)
            if url.scheme or url.netloc:
                continue
            target = ((path.parent / unquote(url.path)) if url.path else path).resolve()
            if not target.is_file():
                errors.append(f"site_en/{rel}: broken link {href}")
                continue
            if target.suffix.lower() == ".html" and url.fragment:
                if target in parsed_pages:
                    target_ids = parsed_pages[target].ids
                else:
                    extra = quality_check.PageParser()
                    extra.feed(target.read_text("utf-8"))
                    target_ids = extra.ids
                if unquote(url.fragment) not in target_ids:
                    errors.append(f"site_en/{rel}: missing anchor {href}")


def check_pdfs(errors: list[str]) -> None:
    max_bytes = 50 * 1024 * 1024  # 50 MiB GitHub warning threshold
    for pdf_name, end_marker in (
        ("microphone-array-tutorial.pdf", "全书完"),
        ("microphone-array-tutorial-en.pdf", "End of Monograph"),
    ):
        pdf_path = DIST / pdf_name
        if not pdf_path.is_file():
            errors.append(f"Missing dist/{pdf_name}")
            continue
        sz = pdf_path.stat().st_size
        if sz >= max_bytes:
            errors.append(f"dist/{pdf_name} exceeds 50 MiB limit: {sz / (1024 * 1024):.2f} MiB")
        reader = PdfReader(str(pdf_path))
        if len(reader.pages) < 500:
            errors.append(f"dist/{pdf_name} page count suspiciously low: {len(reader.pages)}")
        build_pdf.validate_pdf_structure_tags(reader)
        build_pdf.validate_pdf_links(pdf_path)
        last_txt = reader.pages[-1].extract_text() or ""
        if end_marker not in last_txt:
            errors.append(f"dist/{pdf_name} final page missing marker {end_marker!r}")


def check_language_bars(errors: list[str]) -> None:
    """Every Chinese/English Markdown pair must carry a current GitHub-visible switch line."""
    from scripts import language_bars

    errors.extend(language_bars.sync(check=True))


def main() -> None:
    errors: list[str] = []
    check_licenses(errors)
    check_language_bars(errors)
    check_paired_corpus(CHAPTERS_ZH, CHAPTERS_EN, "chapters_en", errors)
    check_paired_corpus(RESEARCH_ZH, RESEARCH_EN, "research_en", errors)
    check_en_repo_docs(errors)
    check_figures_and_svgs(errors)
    check_site_en(errors)
    check_pdfs(errors)
    if errors:
        print(f"ENGLISH QUALITY CHECK FAILED ({len(errors)} issues):", file=sys.stderr)
        for err in errors[:60]:
            print(f"  - {err}", file=sys.stderr)
        raise SystemExit(1)
    print("ENGLISH QUALITY CHECK PASSED")


if __name__ == "__main__":
    main()
