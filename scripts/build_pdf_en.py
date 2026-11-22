#!/usr/bin/env python3
"""Build the English A4 Monograph PDF (`dist/microphone-array-tutorial-en.pdf`)
from `chapters_en/` (19 files) and `figures_en/` (87 figures).

Usage (from repository root):
    .venv/bin/python scripts/build_pdf_en.py
    .venv/bin/python scripts/build_pdf_en.py --html-only
    .venv/bin/python scripts/build_pdf_en.py --pdf-only
"""
from __future__ import annotations

import argparse
import functools
import hashlib
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time
from xml.sax.saxutils import escape as xml_escape

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import markdown
from pypdf import PdfReader, PdfWriter
from pypdf.generic import (
    ArrayObject,
    DictionaryObject,
    Fit,
    NameObject,
    TextStringObject,
)

from scripts import build_pdf, build_site, build_site_en
from scripts.build_markdown_helpers import shield_math, unshield_math, validate_url_schemes
from scripts.chapter_identity import (
    LEGACY_HTML_ROUTES,
    canonical_exercise,
    canonical_source,
    exercise_alias_records,
    inject_exercise_aliases,
)
from scripts.heading_aliases import has_historical_sequential_aliases, historical_aliases

SRC_EN = ROOT / "chapters_en"
SRC_ZH = ROOT / "chapters"
FIG_EN = ROOT / "figures_en"
OUT = ROOT / "dist"

CHAPTERS_EN = [
    ("00_overview.md", "Course Overview & Reading Guide"),
    *build_site_en.CHAPTERS_EN,
]

HTML_TO_CH = {fname.replace(".md", ".html"): build_pdf.chapter_number(fname) for fname, _ in CHAPTERS_EN}
HTML_TO_CH.update({old: HTML_TO_CH[current] for old, current in LEGACY_HTML_ROUTES.items()})

CSS_EN = (
    build_pdf.CSS.replace(
        'font-family:"STHeiti","Hiragino Sans GB","Microsoft YaHei",sans-serif',
        'font-family:"Helvetica Neue",Arial,-apple-system,"Segoe UI",sans-serif',
    )
    .replace(
        'font-family:Menlo,Monaco,Consolas,"STHeiti","Hiragino Sans GB","Microsoft YaHei",monospace',
        'font-family:Menlo,Monaco,Consolas,monospace',
    )
    .replace("line-height:1.74", "line-height:1.66")
)


def source_digest_en() -> str:
    digest = hashlib.sha256()
    paths = sorted(SRC_EN.glob("*.md"))
    paths += [build_site_en.RESEARCH_EN / name for name, _ in build_site_en.RESEARCH_EN_PAGES]
    paths += [
        Path(__file__),
        ROOT / "scripts" / "build_pdf.py",
        ROOT / "scripts" / "build_site_en.py",
        ROOT / "scripts" / "build_markdown_helpers.py",
        build_site.INLINE_LAYOUT_PATH,
        build_pdf.PDF_MATH_SCRIPT,
    ]
    paths += sorted(FIG_EN.glob("fig*.png"))
    for path in paths:
        digest.update(path.relative_to(ROOT).as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()[:12]


def is_exercise_anchor_en(fragment: str, source_path: Path) -> bool:
    if not re.fullmatch(r"e\d{2}-\d{2}", fragment):
        return False
    name = canonical_source(source_path.name)
    canon = canonical_exercise(name, fragment)
    en_file = SRC_EN / name
    return en_file.is_file() and bool(
        re.search(rf'<a\s+id=[\"\']{re.escape(canon)}[\"\']\s*>', en_file.read_text(encoding="utf-8"))
    )


def is_explicit_anchor_en(fragment: str, source_path: Path) -> bool:
    if not re.fullmatch(r"[a-z0-9][a-z0-9_-]*", fragment):
        return False
    name = canonical_source(source_path.name)
    en_file = SRC_EN / name
    return en_file.is_file() and bool(
        re.search(rf'<a\s+id=[\"\']{re.escape(fragment)}[\"\']\s*>', en_file.read_text(encoding="utf-8"))
    )


@functools.lru_cache(maxsize=None)
def _chapter_fragment_map(name: str) -> dict[str, str]:
    en_path = SRC_EN / name
    zh_path = SRC_ZH / name
    if not en_path.is_file() or not zh_path.is_file():
        return {}
    en_records = build_site.heading_records(build_site.parse_headings(en_path.read_text(encoding="utf-8")))
    zh_records = build_site.heading_records(build_site.parse_headings(zh_path.read_text(encoding="utf-8")))
    mapping: dict[str, str] = {}
    for idx, (_lvl, _txt, en_primary, en_legacy) in enumerate(en_records):
        target = "sec-1" if idx == 0 else en_primary
        mapping[en_primary] = target
        if not has_historical_sequential_aliases(name):
            mapping.setdefault(en_legacy, target)
        for old in historical_aliases(name, en_primary):
            mapping.setdefault(old, target)
        if idx < len(zh_records):
            zh_primary = zh_records[idx][2]
            mapping.setdefault(zh_primary, target)
            for old in historical_aliases(name, zh_primary):
                mapping.setdefault(old, target)
    return mapping


def _canonical_en_fragment(md_name: str, fragment: str) -> str:
    """Resolve paired Chinese or historical heading fragments to the exact English primary heading ID."""
    mapping = _chapter_fragment_map(canonical_source(md_name))
    return mapping.get(fragment, fragment)


def rewrite_book_links_en(html: str) -> str:
    def repl(m):
        fname, fragment, text = m.group(1), m.group(2), m.group(3)
        idx = HTML_TO_CH.get(fname)
        if idx is None:
            return m.group(0)
        md_name = "00_overview.md" if fname == "index.html" else fname.replace(".html", ".md")
        if fragment:
            fragment = _canonical_en_fragment(md_name, fragment)
        if not fragment or fragment == "sec-1":
            target = f"ch-{idx}"
        elif idx == 0 and fragment == "course-cases":
            target = "ch-0-audio-assets"
        elif re.fullmatch(r"sec-(?:\d+(?:-\d+)*|u-[0-9a-f]{10}(?:-\d+)?)", fragment):
            target = f"ch-{idx}-{fragment}"
        elif is_exercise_anchor_en(fragment, SRC_EN / md_name):
            canon = canonical_exercise(canonical_source(md_name), fragment)
            target = f"ch-{idx}-{canon}"
        elif is_explicit_anchor_en(fragment, SRC_EN / md_name):
            target = f"ch-{idx}-{fragment}"
        else:
            raise ValueError(f"English PDF cannot map chapter anchor: {fname}#{fragment}")
        return f'<a href="#{target}">{text}</a>'

    return re.sub(r'<a href="(?:\./)?([^"#/]+\.html)(?:#([^"]+))?">(.*?)</a>', repl, html, flags=re.S)


def rewrite_repository_links_en(html: str, source_path: Path) -> str:
    chapters_en_map = {(SRC_EN / name).resolve(): build_pdf.chapter_number(name) for name, _ in CHAPTERS_EN}
    chapters_zh_map = {(SRC_ZH / name).resolve(): build_pdf.chapter_number(name) for name, _ in CHAPTERS_EN}

    def transform(href: str) -> str:
        parsed, target = build_site.local_link_target(href, source_path)
        if target is not None and target.parent in {SRC_EN.resolve(), SRC_ZH.resolve()}:
            target = (SRC_EN / canonical_source(target.name)).resolve()
        if not parsed.scheme and not parsed.netloc and not parsed.path and not parsed.query and parsed.fragment:
            target = Path(source_path).resolve()
        if target is None:
            return href
        idx = chapters_en_map.get(target) if target in chapters_en_map else chapters_zh_map.get(target)
        if idx is not None and not parsed.query:
            fragment = _canonical_en_fragment(target.name, parsed.fragment) if parsed.fragment else ""
            if not fragment or fragment == "sec-1":
                return f"#ch-{idx}"
            if idx == 0 and fragment == "course-cases":
                return "#ch-0-audio-assets"
            if re.fullmatch(r"sec-(?:\d+(?:-\d+)*|u-[0-9a-f]{10}(?:-\d+)?)", fragment):
                return f"#ch-{idx}-{fragment}"
            if is_exercise_anchor_en(fragment, SRC_EN / target.name):
                canon = canonical_exercise(target.name, fragment)
                return f"#ch-{idx}-{canon}"
            if is_explicit_anchor_en(fragment, SRC_EN / target.name):
                return f"#ch-{idx}-{fragment}"
            raise ValueError(f"English PDF cannot map chapter anchor: {target.name}#{fragment}")
        return build_site.repository_url(parsed, target)

    return build_site.rewrite_href_targets(html, transform)


def remove_page_info_en(html: str) -> str:
    html = re.sub(r"<blockquote>\s*<p>📄.*?</p>\s*</blockquote>", "", html, flags=re.S)
    return re.sub(r"<p>📄\s*(?:Page Info|Document Info|本篇信息).*?</p>", "", html, flags=re.S)


def append_book_end_en(html: str) -> str:
    pattern = r"(<p>(?:(?!<p>).)*?</p>\s*(?:<hr\s*/?>\s*)?)$"
    result, count = re.subn(
        pattern,
        r'<div class="book-ending">\1<div class="book-end">End of Monograph</div></div>',
        html,
        flags=re.S,
    )
    if count != 1:
        raise ValueError("Final chapter must end with a normal paragraph <p>...</p> for End of Monograph")
    return result


def build_html_en(build_date: str | None = None):
    build_pdf.check_mathjax_assets()
    body_parts = []
    outline = []
    n_imgs = 0
    for position, (fname, label) in enumerate(CHAPTERS_EN):
        i = build_pdf.chapter_number(fname)
        en_md = (SRC_EN / fname).read_text(encoding="utf-8")
        zh_md = (SRC_ZH / fname).read_text(encoding="utf-8")
        en_records = build_site.heading_records(build_site.parse_headings(en_md))
        zh_records = build_site.heading_records(build_site.parse_headings(zh_md))
        md, repo = shield_math(en_md)
        html = markdown.markdown(md, extensions=["tables", "fenced_code", "sane_lists"])
        html = unshield_math(html, repo)
        validate_url_schemes(html)

        source_heading = [0]
        emitted_chapter_ids: set[str] = set()

        def tag_source_heading(m):
            h_idx = source_heading[0]
            source_heading[0] += 1
            primary = en_records[h_idx][2] if h_idx < len(en_records) else build_pdf.heading_anchor(m.group(2), source_heading[0])
            emitted_chapter_ids.add(primary)
            legacy = f"sec-{source_heading[0]}"
            alias_candidates = []
            html_anchor = build_pdf.heading_anchor(m.group(2), source_heading[0])
            if html_anchor != primary:
                alias_candidates.append(html_anchor)
            if primary != legacy and not has_historical_sequential_aliases(fname):
                alias_candidates.append(legacy)
            for old in historical_aliases(fname, primary):
                if old != primary:
                    alias_candidates.append(old)
            if h_idx < len(zh_records):
                zh_primary = zh_records[h_idx][2]
                if zh_primary != primary:
                    alias_candidates.append(zh_primary)
                for old in historical_aliases(fname, zh_primary):
                    if old != primary:
                        alias_candidates.append(old)
            alias_spans = []
            for cand in alias_candidates:
                if cand not in emitted_chapter_ids:
                    emitted_chapter_ids.add(cand)
                    alias_spans.append(f'<span id="ch-{i}-{cand}" class="anchor-alias"></span>')
            return f'{"".join(alias_spans)}<{m.group(1)} id="ch-{i}-{primary}">{m.group(2)}</{m.group(1)}>'

        html = re.sub(
            r'(<a\s+id=[\"\'])([a-z0-9][a-z0-9_-]*)([\"\']\s*>)',
            lambda m: f"{m.group(1)}ch-{i}-{m.group(2)}{m.group(3)}",
            html,
        )
        html = inject_exercise_aliases(html, fname, prefix=f"ch-{i}-")
        html = re.sub(r"<(h[1-4])>(.*?)</\1>", tag_source_heading, html, flags=re.S)
        html = rewrite_repository_links_en(html, SRC_EN / fname)
        html = rewrite_book_links_en(html)
        html = re.sub(r"^\s*(?:<blockquote>.*?</blockquote>\s*)+(?:<hr\s*/?>\s*)?", "", html, flags=re.S)
        html = remove_page_info_en(html)
        html = build_pdf.strip_terminal_rule(html)
        html = re.sub(
            r'<h[12][^>]*\bid="([^"]+)"[^>]*>.*?</h[12]>',
            r'<span id="\1" class="anchor-alias"></span>',
            html,
            count=1,
            flags=re.S,
        )
        if fname != "00_overview.md":
            html = re.sub(
                r"<h([34])([^>]*)>(.*?)</h\1>",
                lambda m: f"<h{int(m.group(1)) - 1}{m.group(2)}>{m.group(3)}</h{int(m.group(1)) - 1}>",
                html,
                flags=re.S,
            )

        secs = []

        def tag_outline_heading(m):
            level, attrs, inner = int(m.group(1)), m.group(2), m.group(3)
            match = re.search(r'\bid="([^"]+)"', attrs)
            if level == 2:
                fallback = f"ch-{i}-s{len(secs)}"
            else:
                parent_index = max(0, len(secs) - 1)
                child_index = len(secs[-1][2]) if secs else 0
                fallback = f"ch-{i}-s{parent_index}-{child_index}"
            section_id = match.group(1) if match else fallback
            if level == 2:
                secs.append((build_pdf.plain_text(inner), section_id, []))
            elif fname in build_pdf.PDF_THIRD_LEVEL_FILES and secs:
                secs[-1][2].append((build_pdf.plain_text(inner), section_id))
            if match:
                return f"<h{level}{attrs}>{inner}</h{level}>"
            return f'<h{level} id="{section_id}"{attrs}>{inner}</h{level}>'

        html = re.sub(r"<h([23])([^>]*)>(.*?)</h\1>", tag_outline_heading, html, flags=re.S)
        # Ensure figures point to ../figures_en/
        html = html.replace('src="../figures/', 'src="../figures_en/')
        n_imgs += len(re.findall(r"<img ", html))
        if position == len(CHAPTERS_EN) - 1:
            html = append_book_end_en(html)
        body_parts.append(f'<div class="chap" id="ch-{i}"><h1>{label}</h1>{html}</div>')
        outline.append((label, f"ch-{i}", secs))

    toc = ['<section class="toc-page"><h1>Contents</h1><ul class="toc">']
    for label, cid, secs in outline:
        toc.append(f'<li><a href="#{cid}">{label}</a>')
        if secs:
            toc.append('<ul class="sec">')
            for title, sid, subsecs in secs:
                toc.append(f'<li><a href="#{sid}">{build_pdf.toc_label(title)}</a>')
                if subsecs:
                    toc.append('<ul class="subsec">')
                    toc.extend(
                        f'<li><a href="#{subid}">{build_pdf.toc_label(subtitle)}</a></li>'
                        for subtitle, subid in subsecs
                    )
                    toc.append("</ul>")
                toc.append("</li>")
            toc.append("</ul>")
        toc.append("</li>")
    toc.append("</ul></section>")

    date_s = build_pdf.resolve_build_date(build_date)
    cover = (
        '<div class="cover">'
        '<article class="cover-banner">'
        '<p class="cover-kicker">OPEN TECHNICAL MONOGRAPH · THEORY, PURE NUMPY &amp; INDUSTRIAL PRACTICE</p>'
        "<h1>Microphone Array Signal Processing</h1>"
        '<p class="sub">From Mathematical Foundations and Core Algorithms to Real-Time Hardware Selection &amp; Frontier Sound Fields (A4 Monograph)</p>'
        "</article>"
        '<table class="cover-table">'
        "<thead><tr>"
        '<th scope="col">Foundations &amp; Six Core Algorithms (Overview, Ch 0–9)</th>'
        '<th scope="col">Real-Time Deployment, Frontier Arrays &amp; Appendices (Ch 10–17)</th>'
        "</tr></thead>"
        "<tbody>"
        "<tr>"
        "<td><strong>Part I · Mathematical Base &amp; Acoustic Modeling (Overview, Ch 0–3)</strong>"
        "SO(3) spherical geometry, complex phasor STFT/COLA, subspace projections, and Wirtinger calculus; microsecond wavefront TDOAs, steering vectors, room acoustics, and linear/circular/spherical/nested coarray limits.</td>"
        "<td><strong>Part III · Real-Time Systems &amp; Pareto Selection (Ch 10–11)</strong>"
        "Front-end cascade topology and AEC loopback taps, crystal sampling-rate offset (SRO), channel calibration and dead-mic fallback, fixed-point MACs/power budgets, and testable acceptance contracts.</td>"
        "</tr>"
        "<tr>"
        "<td><strong>Part II · Six Core Array Algorithms (Ch 4–9)</strong>"
        "DOA localization, fixed/adaptive/neural beamforming, full-duplex acoustic echo cancellation (AEC), WPE blind dereverberation, multichannel speech separation, and Bayesian multi-target tracking.</td>"
        "<td><strong>Part IV · Frontier Chapters &amp; Reference Appendices (Ch 12–17)</strong>"
        "Deconvolution acoustic camera imaging, low-bandwidth wireless distributed arrays, Ambisonics &amp; wearable micro-arrays, conditioned target extraction &amp; diarization, plus master notation lookup and triage trees.</td>"
        "</tr>"
        "</tbody></table>"
        '<table class="cover-stats-table">'
        "<thead><tr>"
        f'<th scope="col">{len(CHAPTERS_EN)} Chapters &amp; Guides</th>'
        '<th scope="col">352 Runnable Exercises</th>'
        '<th scope="col">87 Traced Figures</th>'
        '<th scope="col">283 Companion WAVs</th>'
        "</tr></thead>"
        "<tbody><tr>"
        "<td>Overview + 16 Chapters + 2 Appendices</td>"
        "<td>Exact Hand Anchors + Pure NumPy</td>"
        "<td>Script-Generated &amp; Numerically Verified</td>"
        "<td>Real Room Recordings &amp; Controlled Fixtures</td>"
        "</tr></tbody></table>"
        f'<article class="meta"><p style="margin:0">Build date {date_s} · Source sha256 {source_digest_en()}<br/>'
        f"{len(CHAPTERS_EN)} parts in one volume · Full-text searchable with 3-level bookmarks and accessible SRE MathSpeak formula tags</p></article>"
        "</div>"
    )

    real_ids = re.findall(r'<(?:h[1-4]|div|a)\b[^>]*\bid="([^"]+)"', "\n".join(body_parts))
    references = (
        '<div class="pdf-anchor-references" aria-hidden="true">'
        + "".join(f'<a href="#{identity}">.</a>' for identity in real_ids)
        + "</div>"
    )
    page = (
        '<!DOCTYPE html><html lang="en" data-tutorial-print-layout="a4"><head><meta charset="utf-8">'
        f"<title>Microphone Array Signal Processing (Full Monograph)</title><style>{CSS_EN}</style>"
        "<script>\n"
        + build_site.INLINE_LAYOUT_PATH.read_text(encoding="utf-8")
        + "\n</script><script>"
        + build_pdf.PDF_MATH_SCRIPT.read_text(encoding="utf-8")
        + "</script>"
        '<script defer src="../scripts/vendor/mathjax-3.2.2/tex-svg.js"></script>'
        "</head><body>"
        + references
        + cover
        + "\n".join(toc)
        + "\n".join(body_parts)
        + "</body></html>"
    )
    page = build_pdf.label_pdf_math_sources(page)
    n_secs = sum(len(s) for _, _, s in outline)
    n_subsecs = sum(len(subsecs) for _, _, secs in outline for _t, _sid, subsecs in secs)
    unique_imgs = len(set(re.findall(r'<img [^>]*src="([^"]+)"', page)))
    print(f"Combined EN: {len(CHAPTERS_EN)} parts / {n_secs} sections / {n_subsecs} subsections / {unique_imgs} unique figures ({n_imgs} refs)")
    return page, outline


def pdf_alias_targets_en() -> dict[str, str]:
    aliases: dict[str, str] = {}

    def bind(alias: str, target: str) -> None:
        if alias == target:
            return
        if alias not in aliases:
            aliases[alias] = target

    for name, _label in CHAPTERS_EN:
        prefix = f"ch-{build_pdf.chapter_number(name)}-"
        en_records = build_site.heading_records(
            build_site.parse_headings((SRC_EN / name).read_text(encoding="utf-8"))
        )
        zh_records = build_site.heading_records(
            build_site.parse_headings((SRC_ZH / name).read_text(encoding="utf-8"))
        )
        for index, (_level, _title, primary, legacy) in enumerate(en_records):
            target = prefix + primary if index else prefix[:-1]
            bind(prefix + primary, target)
            if not has_historical_sequential_aliases(name):
                bind(prefix + legacy, target)
            for old in historical_aliases(name, primary):
                bind(prefix + old, target)
            if index < len(zh_records):
                zh_primary = zh_records[index][2]
                bind(prefix + zh_primary, target)
                for old in historical_aliases(name, zh_primary):
                    bind(prefix + old, target)
        for current, record in exercise_alias_records(name).items():
            resolved_target = aliases.get(prefix + current, prefix + current)
            for old_id in record.get("legacy_ids", [record["original_id"]]):
                bind(prefix + old_id, resolved_target)
    return aliases


def _verified_contourless_space_font(font_ref) -> bool:
    from io import BytesIO
    from fontTools.ttLib import TTFont, TTLibError

    if build_pdf.verified_layout_space_font(font_ref):
        return True
    try:
        font = font_ref.get_object()
        if font.get("/Subtype") != "/Type0" or font.get("/Encoding") != "/Identity-H":
            return False
        descendant = font["/DescendantFonts"][0].get_object()
        if descendant.get("/Subtype") != "/CIDFontType2" or descendant.get("/CIDToGIDMap") != "/Identity":
            return False
        data = descendant["/FontDescriptor"]["/FontFile2"].get_data()
        with TTFont(BytesIO(data)) as embedded:
            glyph = embedded.getGlyphOrder()[3]
            return "glyf" in embedded and embedded["glyf"][glyph].numberOfContours == 0
    except (KeyError, IndexError, ValueError, TypeError, OSError, TTLibError):
        return False


def _mark_orphan_layout_span_spaces(writer: PdfWriter) -> int:
    """Wrap top-level `/Span << /ActualText ( ) >> BDC ... (\\000\\003) Tj EMC` blocks
    that lie outside any `/MCID` or `/Artifact` and use a verified contourless space font
    as `/Artifact << /Type /Layout >> BDC ... EMC`.
    """
    import zlib
    from pypdf.generic import EncodedStreamObject, IndirectObject, StreamObject

    token_re = re.compile(rb"(/[A-Za-z0-9_]+)[ \t\r\n]*<<([^>]*)>>[ \t\r\n]*BDC|\b(EMC)\b")
    span_space_re = re.compile(
        rb"^/Span\s*<<\s*/ActualText\s*\(\s*\)\s*>>\s*BDC\s*"
        rb"(?:/(F\d+)\s+[0-9.]+\s+Tf\s+)?"
        rb"(?:[0-9.+-]+\s+){2,6}(?:Tm|Td)\s+"
        rb"\(\\000\\003\)\s*Tj\s*EMC$",
        re.S,
    )
    font_proofs: dict[tuple, bool] = {}
    fixed_total = 0

    for page in writer.pages:
        ref = page.raw_get("/Contents")
        if not hasattr(ref, "idnum"):
            continue
        cid = ref.idnum
        obj = writer._objects[cid - 1]
        if not isinstance(obj, StreamObject) or obj.get("/DecodeParms") is not None:
            continue
        raw = obj.get_data()
        if b"/Span" not in raw:
            continue

        # Locate top-level /Span BDC..EMC ranges outside any /MCID or /Artifact
        stack = []
        replacements = []
        for m in token_re.finditer(raw):
            if m.group(1):
                tag, props = m.group(1), m.group(2)
                stack.append((tag, b"/MCID" in props, tag == b"/Artifact", m.start()))
            else:
                if not stack:
                    continue
                tag, has_mcid, is_art, start_pos = stack.pop()
                if tag == b"/Span" and not has_mcid and not is_art and not any(h or a for _, h, a, _ in stack):
                    end_pos = m.end()
                    block = raw[start_pos:end_pos]
                    sm = span_space_re.match(block)
                    if sm:
                        fname = "/" + sm.group(1).decode("ascii") if sm.group(1) else None
                        ok = True
                        if fname:
                            res = page.get("/Resources", {}).get("/Font", {}).get(fname)
                            if res is not None:
                                resolved = res.get_object()
                                iref = getattr(resolved, "indirect_reference", None)
                                key = ("ind", iref.idnum) if iref else ("dir", id(resolved))
                                if key not in font_proofs:
                                    font_proofs[key] = _verified_contourless_space_font(res)
                                ok = font_proofs[key]
                        if ok:
                            new_block = re.sub(
                                rb"^/Span\s*<<\s*/ActualText\s*\(\s*\)\s*>>\s*BDC",
                                b"/Artifact << /Type /Layout >> BDC",
                                block,
                                count=1,
                            )
                            replacements.append((start_pos, end_pos, new_block))

        if replacements:
            out_parts = []
            pos = 0
            for s, e, nb in replacements:
                out_parts.append(raw[pos:s])
                out_parts.append(nb)
                pos = e
            out_parts.append(raw[pos:])
            new_raw = b"".join(out_parts)
            comp = zlib.compress(new_raw, 9)
            enc = EncodedStreamObject()
            for k, v in obj.items():
                if k not in ("/Length", "/Filter", "/DecodeParms"):
                    enc[NameObject(k)] = v
            enc._data = comp
            enc[NameObject("/Filter")] = NameObject("/FlateDecode")
            enc[NameObject("/Length")] = build_pdf.pypdf_number(len(comp)) if hasattr(build_pdf, "pypdf_number") else __import__("pypdf.generic", fromlist=["NumberObject"]).NumberObject(len(comp))
            enc.indirect_reference = obj.indirect_reference or IndirectObject(cid, 0, writer)
            writer._objects[cid - 1] = enc
            fixed_total += len(replacements)
    return fixed_total


def normalize_pdf_semantics_en(writer: PdfWriter, outline, aliases=None) -> None:
    _mark_orphan_layout_span_spaces(writer)
    tree = writer.root_object.get("/StructTreeRoot")
    if tree is None:
        raise SystemExit("English PDF missing /StructTreeRoot")
    tree = tree.get_object()
    role_map = tree.get("/RoleMap", DictionaryObject()).get_object()
    for role in ("/Strong", "/Em"):
        role_map[NameObject(role)] = NameObject("/Span")
    tree[NameObject("/RoleMap")] = role_map

    def visit(value):
        node = value.get_object() if hasattr(value, "get_object") else value
        if isinstance(node, list):
            for child in node:
                visit(child)
            return
        if not isinstance(node, dict):
            return
        if node.get("/S") == "/LI":
            children = node.get("/K", [])
            children = children if isinstance(children, list) else [children]
            resolved = [item.get_object() if hasattr(item, "get_object") else item for item in children]
            types = [item.get("/S") if isinstance(item, dict) else None for item in resolved]
            if any(kind not in ("/Lbl", "/LBody") for kind in types):
                body_children = [item for item, kind in zip(children, types) if kind != "/Lbl"]
                body = DictionaryObject({
                    NameObject("/Type"): NameObject("/StructElem"),
                    NameObject("/S"): NameObject("/LBody"),
                    NameObject("/P"): node.indirect_reference,
                    NameObject("/K"): ArrayObject(body_children),
                })
                if "/Pg" in node:
                    body[NameObject("/Pg")] = node.raw_get("/Pg")
                body_ref = writer._add_object(body)
                for item in body_children:
                    item.get_object()[NameObject("/P")] = body_ref
                labels = [item for item, kind in zip(children, types) if kind == "/Lbl"]
                node[NameObject("/K")] = ArrayObject(labels + [body_ref])
        if "/K" in node:
            visit(node["/K"])

    visit(tree.get("/K", []))
    real_targets = {str(key).removeprefix("/") for key in writer.named_destinations}
    titles = {"ch-0-audio-assets": "Audio Suites and Interactive Cases"}
    for chapter_title, chapter_id, sections in outline:
        titles[chapter_id] = chapter_title
        for title, identity, subsections in sections:
            titles[identity] = title
            titles.update({sub_id: sub_title for sub_title, sub_id in subsections})
    for alias, target in (aliases or {}).items():
        if target in titles:
            titles[alias] = titles[target]

    for page in writer.pages:
        page[NameObject("/Tabs")] = NameObject("/S")
        for reference in page.get("/Annots", []):
            annotation = reference.get_object()
            if annotation.get("/Subtype") != "/Link" or "/StructParent" not in annotation or annotation.get("/Contents"):
                continue
            action = annotation.get("/A", {})
            if action.get("/S") == "/URI" and action.get("/URI"):
                description = "Open external link: " + str(action["/URI"])
            else:
                target = str(annotation.get("/Dest", "")).removeprefix("/")
                if target not in real_targets:
                    raise SystemExit(f"English PDF internal link target missing: {target}")
                if target in titles:
                    description = "Jump to: " + titles[target]
                else:
                    exercise = re.search(r"E\d{2}-\d{2}$", (aliases or {}).get(target, target), re.I)
                    description = ("View exercise: " + exercise.group().upper()) if exercise else f"Jump to section: {target}"
            annotation[NameObject("/Contents")] = TextStringObject(description)

    writer.pdf_header = "%PDF-1.7"
    writer.xmp_metadata = (
        '<?xpacket begin="\ufeff" id="W5M0MpCehiHzreSzNTczkc9d"?>'
        '<x:xmpmeta xmlns:x="adobe:ns:meta/"><rdf:RDF '
        'xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">'
        '<rdf:Description rdf:about="" xmlns:dc="http://purl.org/dc/elements/1.1/" '
        'xmlns:xmp="http://ns.adobe.com/xap/1.0/" '
        'xmp:CreatorTool="scripts/build_pdf_en.py">'
        '<dc:title><rdf:Alt><rdf:li xml:lang="x-default">'
        + xml_escape(str(writer.metadata.get("/Title", "Microphone Array Signal Processing")))
        + "</rdf:li></rdf:Alt></dc:title>"
        "<dc:language><rdf:Bag><rdf:li>en-US</rdf:li></rdf:Bag></dc:language>"
        "</rdf:Description></rdf:RDF></x:xmpmeta><?xpacket end=\"w\"?>"
    ).encode("utf-8")
    metadata = writer.root_object["/Metadata"]
    metadata[NameObject("/Type")] = NameObject("/Metadata")
    metadata[NameObject("/Subtype")] = NameObject("/XML")


def print_pdf_en(combined: Path, pdf: Path, timeout_min_pages: int = 100) -> int:
    chrome = build_pdf.find_chrome()
    if not Path(chrome).exists():
        raise SystemExit(f"Chrome not found: {chrome}")
    fd, tmp_name = tempfile.mkstemp(prefix="printed-en-", suffix=".pdf", dir=pdf.parent)
    os.close(fd)
    tmp_pdf = Path(tmp_name)
    tmp_pdf.unlink()
    completed_by_stability = False
    log_text = ""
    try:
        with tempfile.TemporaryDirectory() as udd, tempfile.TemporaryFile(mode="w+") as log:
            proc = subprocess.Popen(
                [
                    chrome,
                    "--headless",
                    "--disable-gpu",
                    "--no-sandbox",
                    "--allow-file-access-from-files",
                    f"--user-data-dir={udd}",
                    "--timeout=240000",
                    "--virtual-time-budget=35000",
                    "--export-tagged-pdf",
                    f"--print-to-pdf={tmp_pdf}",
                    "--no-pdf-header-footer",
                    combined.as_uri(),
                ],
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
            )
            deadline, last_size, stable_since = time.time() + 660, -1, None
            while time.time() < deadline:
                if proc.poll() is not None:
                    break
                if tmp_pdf.exists():
                    sz = tmp_pdf.stat().st_size
                    now = time.time()
                    if sz == last_size and sz > 64 * 1024:
                        if stable_since is not None and now - stable_since > 15:
                            completed_by_stability = True
                            break
                    else:
                        last_size, stable_since = sz, now
                time.sleep(5)
            if proc.poll() is None:
                proc.kill()
                proc.wait()
            log.seek(0)
            log_text = log.read()[-1000:]
            if proc.returncode not in (0, None) and not completed_by_stability:
                raise SystemExit(f"Chrome print failed ({proc.returncode})\n{log_text}")
        if not tmp_pdf.exists():
            raise SystemExit(f"English PDF not generated\n{log_text}")
        print("generated raw EN", tmp_pdf, tmp_pdf.stat().st_size // 1024 // 1024, "MiB")
        writer = PdfWriter(clone_from=PdfReader(str(tmp_pdf)))
        build_pdf.tag_pdf_formulas(writer, build_pdf.pdf_math_sources(combined.read_text(encoding="utf-8")))
        with tmp_pdf.open("wb") as stream:
            writer.write(stream)
        reader = PdfReader(str(tmp_pdf))
        build_pdf.validate_pdf_structure_tags(reader)
        npages = len(reader.pages)
        page_texts = []
        scale_samples = []
        for page_number, page in enumerate(reader.pages, 1):
            def _visit_scale(text, cm, tm, _font, font_size, _pno=page_number):
                if (
                    abs(font_size - 16) < 0.01
                    and len(re.findall(r"[A-Za-z]", text)) >= 12
                    and abs(cm[1]) < 1e-6
                    and abs(tm[1]) < 1e-6
                ):
                    scale_samples.append((_pno, abs(cm[0] * tm[0])))
            page_texts.append(page.extract_text(visitor_text=_visit_scale) or "")
        if npages < timeout_min_pages:
            raise SystemExit(f"English PDF page count too small ({npages})")
        if "End of Monograph" not in page_texts[-1]:
            raise SystemExit("English PDF final page missing 'End of Monograph'")
        raw_math_pages = [
            (i + 1, page_text)
            for i, page_text in enumerate(page_texts)
            if build_pdf.contains_unrendered_math(page_text)
        ]
        if raw_math_pages:
            samples = [
                f"p{number}: "
                + " | ".join(line for line in page_text.splitlines() if build_pdf.contains_unrendered_math(line))
                for number, page_text in raw_math_pages[:5]
            ]
            raise SystemExit("English PDF text layer has unrendered TeX:\n" + "\n".join(samples))
        build_pdf.validate_pdf_text_codepoints("\n".join(page_texts))
        if not scale_samples:
            raise SystemExit("English PDF missing 16px body scale samples")
        too_small = [(p, round(s, 4)) for p, s in scale_samples if s < 0.74]
        if too_small:
            raise SystemExit(f"English PDF body scaled down (<0.75): {too_small[:5]}")
        os.replace(tmp_pdf, pdf)
    finally:
        tmp_pdf.unlink(missing_ok=True)
    return npages


def add_bookmarks_en(pdf_path: Path, outline, aliases=None):
    reader = PdfReader(str(pdf_path))
    writer = PdfWriter(clone_from=reader)
    writer.add_metadata({
        "/Title": "Microphone Array Signal Processing",
        "/Subject": "From Fundamentals to Industrial Practice",
        "/Creator": "scripts/build_pdf_en.py",
        "/SourceDigest": source_digest_en(),
    })
    writer._root_object.update({NameObject("/Lang"): TextStringObject("en-US")})
    if aliases is not None:
        build_pdf.add_pdf_named_aliases(reader, writer, aliases)
    normalize_pdf_semantics_en(writer, outline, aliases)
    n_ch, n_sec, n_subsec = 0, 0, 0

    def add_heading(title, heading_id, parent=None):
        destination = build_pdf.heading_destination(reader, heading_id)
        target = destination.dest_array
        return writer.add_outline_item(
            title,
            reader.get_destination_page_number(destination),
            parent=parent,
            fit=Fit(str(target[1]), tuple(target[2:])),
        )

    for label, cid, secs in outline:
        parent = add_heading(label, cid)
        n_ch += 1
        for title, sid, subsecs in secs:
            section_parent = add_heading(title, sid, parent)
            n_sec += 1
            for subtitle, subid in subsecs:
                add_heading(subtitle, subid, section_parent)
                n_subsec += 1

    fd, tmp_name = tempfile.mkstemp(prefix="bookmarked-en-", suffix=".pdf", dir=pdf_path.parent)
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        with open(tmp_path, "wb") as f:
            build_pdf.write_pdf_with_object_streams(writer, f)
        tagged = PdfReader(str(tmp_path))
        build_pdf.validate_pdf_structure_tags(tagged)
        if aliases is not None:
            build_pdf.validate_pdf_named_aliases(tagged, aliases)
        os.replace(tmp_path, pdf_path)
    finally:
        tmp_path.unlink(missing_ok=True)
    print(f"Bookmarks EN: {n_ch} parts, {n_sec} sections, {n_subsec} subsections")
    return n_ch, n_sec, n_subsec


def check_figures_en() -> None:
    missing = []
    refs = set()
    for fname, _ in CHAPTERS_EN:
        md = (SRC_EN / fname).read_text(encoding="utf-8")
        for m in re.finditer(r"\.\./figures(?:_en)?/([^\")\s]+)", md):
            refs.add(m.group(1))
            path = FIG_EN / m.group(1)
            if not path.exists() or path.stat().st_size == 0:
                missing.append(f"{fname}: {m.group(1)}")
    if missing:
        raise SystemExit("Missing English figures:\n" + "\n".join(missing))
    if len(refs) != 87:
        raise SystemExit(f"Expected 87 unique figure references in chapters_en/, got {len(refs)}")
    print(f"English figures check passed ({len(CHAPTERS_EN)} parts, {len(refs)} unique figures)")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="Build English Monograph PDF: chapters_en/ -> dist/microphone-array-tutorial-en.pdf")
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--html-only", action="store_true")
    mode.add_argument("--pdf-only", action="store_true")
    ap.add_argument("--no-bookmarks", action="store_true")
    ap.add_argument("--min-pages", type=int, default=100)
    ap.add_argument("--build-date")
    args = ap.parse_args(argv)

    build_pdf.check_mathjax_assets()
    OUT.mkdir(exist_ok=True)
    combined = OUT / "combined_en.html"
    pdf = OUT / "microphone-array-tutorial-en.pdf"
    outline = None
    tmp_html = None
    tmp_pdf = None
    try:
        if not args.pdf_only:
            check_figures_en()
            page, outline = build_html_en(args.build_date)
            fd, tmp_name = tempfile.mkstemp(prefix="combined-en-", suffix=".html", dir=OUT)
            os.close(fd)
            tmp_html = Path(tmp_name)
            tmp_html.write_text(page, encoding="utf-8")
        if args.html_only:
            build_site.publish_files([(tmp_html, combined)])
            print("saved", combined, "(html only)")
            return
        if args.pdf_only:
            html = combined.read_text(encoding="utf-8")
            outline = build_pdf.outline_from_html(html)
        fd, tmp_name = tempfile.mkstemp(prefix="validated-en-", suffix=".pdf", dir=OUT)
        os.close(fd)
        tmp_pdf = Path(tmp_name)
        print_pdf_en(tmp_html or combined, tmp_pdf, args.min_pages)
        if not args.no_bookmarks:
            add_bookmarks_en(tmp_pdf, outline, pdf_alias_targets_en())
        build_pdf.validate_pdf_links(tmp_pdf)
        replacements = [] if tmp_html is None else [(tmp_html, combined)]
        replacements.append((tmp_pdf, pdf))
        build_site.publish_files(replacements)
        os.chmod(pdf, 0o644)
        sz_mib = pdf.stat().st_size / (1024 * 1024)
        print(f"saved {pdf} ({sz_mib:.2f} MiB)")
    finally:
        if tmp_html is not None:
            tmp_html.unlink(missing_ok=True)
        if tmp_pdf is not None:
            tmp_pdf.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
