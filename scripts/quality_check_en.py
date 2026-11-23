#!/usr/bin/env python3
"""Comprehensive quality gate for the English Edition, Dual License, and <50 MiB PDFs."""
from __future__ import annotations

from collections import Counter
from html import unescape
from html.parser import HTMLParser
import hashlib
import json
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


# Independent published English titles; never derive expected labels from a PDF.
EXPECTED_PDF_CHAPTERS_EN = (('00_overview.md', 'Course Overview & Reading Guide'),
 ('00_math-foundations.md', 'Chapter 0 · Mathematical & Signal Analysis Foundations'),
 ('01_problem-definition.md', 'Chapter 1 · Problem Formulation & Binaural Cues'),
 ('02_basics-signal-model.md', 'Chapter 2 · Acoustic Propagation & Array Signal Model'),
 ('03_array-geometry.md', 'Chapter 3 · Microphone Array Geometry'),
 ('04_doa-estimation.md', 'Chapter 4 · Sound Source Localization (DOA)'),
 ('05_beamforming.md', 'Chapter 5 · Beamforming'),
 ('06_aec.md', 'Chapter 6 · Acoustic Echo Cancellation (AEC)'),
 ('07_wpe-dereverberation.md', 'Chapter 7 · Dereverberation (WPE)'),
 ('08_speech-separation.md', 'Chapter 8 · Multichannel Speech Separation'),
 ('09_source-tracking.md', 'Chapter 9 · Sound Source Tracking'),
 ('10_engineering-practice.md', 'Chapter 10 · Real-Time Systems, Evaluation & Deployment'),
 ('11_selection-guide.md', 'Chapter 11 · System Synthesis & Selection Guide'),
 ('12_acoustic-imaging.md', 'Chapter 12 · Acoustic Imaging & Noise Diagnostics'),
 ('13_distributed-enhancement.md', 'Chapter 13 · Distributed Wireless Microphone Arrays'),
 ('14_spatial-wearable-arrays.md', 'Chapter 14 · Spatial Sound Fields & Wearable Arrays'),
 ('15_target-extraction-diarization.md', 'Chapter 15 · Conditioned Target Extraction & Diarization'),
 ('16_appendix-symbols-math.md', 'Chapter 16 · Appendix A: Symbols, Glossary & Math Lookup'),
 ('17_appendix-guide.md', 'Chapter 17 · Appendix B: Literature Map, Triage & Exercises'))


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


def figure_source_script(filename: str) -> str:
    number = int(re.fullmatch(r"fig(\d{2})_[a-z0-9_]+\.png", filename).group(1))
    singles = {68: "beamforming", 69: "reference", 70: "delay", 71: "css",
               72: "tracking", 73: "channel", 74: "selection"}
    if number in singles:
        return f"scripts/make_{singles[number]}_figures.py"
    if 76 <= number <= 79:
        return "scripts/make_wearable_figures.py"
    if 80 <= number <= 83:
        return "scripts/make_tse_figures.py"
    if 84 <= number <= 87:
        return "scripts/make_math_figures.py"
    if number in (25, 26, 27, 32, 33, 34):
        return "scripts/make_aec_figures.py"
    return "scripts/make_figures.py"


def english_png_issues(path: Path, filename: str) -> list[str]:
    """Fully decode and verify identity without deriving it from the PNG itself."""
    from PIL import Image
    issues = []
    try:
        with Image.open(path) as picture:
            picture.load()
            if picture.format != "PNG":
                raise ValueError("not PNG")
            metadata = dict(picture.info)
            if picture.width < 800 or picture.height < 300:
                issues.append("resolution below publication minimum")
        source = figure_source_script(filename)
        expected = {
            "FigureFilename": filename,
            "SourceScript": source,
            "SourceScriptDigest": hashlib.sha256((ROOT / source).read_bytes()).hexdigest(),
            "TranslationScript": "scripts/make_figures_en.py",
            "TranslationScriptSHA256": hashlib.sha256((ROOT / "scripts/make_figures_en.py").read_bytes()).hexdigest(),
            "TranslationTableSHA256": hashlib.sha256((ROOT / "scripts/figure_en_translations.json").read_bytes()).hexdigest(),
            "TranslationStyleSHA256": hashlib.sha256((ROOT / "scripts/make_figures.py").read_bytes()).hexdigest(),
        }
        for key, value in expected.items():
            if metadata.get(key) != value:
                issues.append(f"{key} mismatch")
        # Translation must retain the source figure's numerical/audio provenance.
        with Image.open(FIGURES_ZH / filename) as original:
            original.load()
            for key in ("GeneratorInputs", "NumericalReportDigest", "AudioManifestDigest"):
                if metadata.get(key) != original.info.get(key):
                    issues.append(f"{key} differs from Chinese numerical figure")
        if "GeneratorInputs" in metadata:
            inputs = json.loads(metadata["GeneratorInputs"])
            if not isinstance(inputs, dict):
                raise ValueError("GeneratorInputs is not an object")
            for relative, digest in inputs.items():
                source_path = ROOT / relative
                if source_path.resolve().is_relative_to(ROOT.resolve()):
                    if hashlib.sha256(source_path.read_bytes()).hexdigest() != digest:
                        issues.append(f"stale generator input: {relative}")
                else:
                    issues.append(f"generator input outside repository: {relative}")
    except (OSError, ValueError, KeyError, TypeError) as exc:
        issues.append(f"invalid PNG or provenance: {exc}")
    return issues


def check_figures_and_svgs(errors: list[str]) -> None:
    from scripts.figure_report_contracts import report_bindings_issues
    for directory in (FIGURES_ZH, FIGURES_EN):
        errors.extend(f"{directory.name}: {issue}" for issue in report_bindings_issues(ROOT, directory))
    expected = set(quality_check.EXPECTED_FIGURES)
    for directory in (FIGURES_ZH, FIGURES_EN):
        actual = {p.name for p in directory.glob("fig*.png")}
        if actual != expected:
            errors.append(f"{directory.name} mismatch: missing={sorted(expected-actual)}, extra={sorted(actual-expected)}")
    for name in sorted(expected):
        errors.extend(f"figures_en/{name}: {issue}"
                      for issue in english_png_issues(FIGURES_EN / name, name))

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


def expected_pdf_outline_en(documents=None):
    """Source heading identities plus the independent published chapter catalog."""
    if documents is None:
        documents = {name: (CHAPTERS_EN / name).read_text(encoding="utf-8")
                     for name, _ in EXPECTED_PDF_CHAPTERS_EN}
    result = []
    for name, label in EXPECTED_PDF_CHAPTERS_EN:
        cid = "ch-" + str(quality_check.EXPECTED_PDF_TOPIC_IDS.get(name, int(name.split("_", 1)[0])))
        level = 2 if name == "00_overview.md" else 3
        sections, current = [], None
        records = quality_check.semantic_heading_ids(documents[name])
        raw = quality_check.markdown_headings(documents[name])
        for (depth, _title, primary), (_depth, title) in zip(records, raw):
            entry = (quality_check.outline_heading_text(title), f"{cid}-{primary}", [])
            if depth == level:
                current = entry
                sections.append(entry)
            elif depth == level + 1 and name in quality_check.EXPECTED_SUBSECTION_COUNTS and current is not None:
                current[2].append(entry[:2])
        if len(sections) != quality_check.EXPECTED_SECTION_COUNTS[name]:
            raise ValueError(f"English source section count differs from publication catalog: {name}")
        if sum(len(item[2]) for item in sections) != quality_check.EXPECTED_SUBSECTION_COUNTS.get(name, 0):
            raise ValueError(f"English source subsection count differs from publication catalog: {name}")
        result.append((label, cid, sections))
    return result


def _outline_identity(outline):
    clean = build_pdf.bookmark_plain_text
    return [(clean(label), cid, [(clean(title), sid, [(clean(sub), subid) for sub, subid in children])
                                for title, sid, children in sections])
            for label, cid, sections in outline]


def cached_html_issues_en(html: str, current_digest: str) -> list[str]:
    digests = re.findall(r"Source sha256 ([0-9a-f]{12})(?![0-9a-f])", html)
    issues = []
    if digests != [current_digest]:
        issues.append("English combined HTML missing, ambiguous or stale source digest")
    try:
        actual = _outline_identity(build_pdf.outline_from_html(html))
        expected = _outline_identity(expected_pdf_outline_en())
        if actual != expected:
            issues.append("English combined HTML outline/heading identities differ from current sources")
        ids = re.findall(r'\bid=["\']([^"\']+)["\']', html)
        if len(ids) != len(set(ids)):
            issues.append("English combined HTML contains duplicate IDs")
    except (KeyError, ValueError, TypeError) as exc:
        issues.append(f"English combined HTML source expectation failed: {exc}")
    return issues


def pdf_outline_issues_en(reader) -> list[str]:
    issues = []
    try:
        expected = expected_pdf_outline_en()
        actual = quality_check.pdf_outline_tree(reader.outline)
        def check_level(nodes, entries, depth=0):
            if len(nodes) != len(entries):
                issues.append(f"English PDF outline length mismatch at depth {depth}")
            for (destination, children), (title, heading_id, subentries) in zip(nodes, entries):
                if destination.title != build_pdf.bookmark_plain_text(title):
                    issues.append(f"English PDF outline title/order mismatch: {heading_id}")
                issues.extend(quality_check.bookmark_destination_issues(reader, destination, heading_id))
                check_level(children, subentries, depth + 1)
        entries = [(label, cid, [(title, sid, [(sub, subid, []) for sub, subid in subs])
                                 for title, sid, subs in sections]) for label, cid, sections in expected]
        check_level(actual, entries)
    except (KeyError, ValueError, TypeError) as exc:
        issues.append(f"English PDF source expectation failed: {exc}")
    return issues


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
        if pdf_name.endswith("-en.pdf"):
            from scripts import build_pdf_en
            digest = build_pdf_en.source_digest_en()
            errors.extend(pdf_outline_issues_en(reader))
            builder = build_pdf_en.build_html_en
        else:
            digest = build_pdf.source_digest()
            builder = build_pdf.build_html
        if reader.metadata is None or reader.metadata.get("/SourceDigest") != digest:
            errors.append(f"dist/{pdf_name} missing or stale SourceDigest")
        try:
            build_pdf.verify_pdf_input_metadata(reader, builder,
                build_pdf_en.source_digest_en if pdf_name.endswith("-en.pdf") else build_pdf.source_digest)
        except (ValueError, TypeError, OSError) as exc:
            errors.append(f"dist/{pdf_name} actual HTML input identity: {exc}")
        if len(reader.pages) < 500:
            errors.append(f"dist/{pdf_name} page count suspiciously low: {len(reader.pages)}")
        build_pdf.validate_pdf_structure_tags(reader)
        build_pdf.validate_pdf_links(pdf_path)
        def _walk_outline(items):
            for it in items:
                if isinstance(it, list):
                    yield from _walk_outline(it)
                else:
                    yield it
        dirty_bookmarks = [
            it.title for it in _walk_outline(reader.outline)
            if "$" in it.title or "\\" in it.title
        ]
        if dirty_bookmarks:
            errors.append(
                f"dist/{pdf_name} contains {len(dirty_bookmarks)} unsanitized TeX bookmarks "
                f"(e.g. {dirty_bookmarks[0]!r})"
            )
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
