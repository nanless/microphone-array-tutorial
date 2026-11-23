#!/usr/bin/env python3
"""Build the English static course site (`site_en/`) from `chapters_en/` and `research_en/`.

Usage (from repository root):
    .venv/bin/python scripts/build_site_en.py

Outputs 33 HTML pages in `site_en/` (`index.html`, 18 chapters/appendices, 8 legacy aliases,
and 6 research pages in `site_en/research/`). Audio links reference `../site/<audio_dir>/*.wav`
so the 283 WAV files are never duplicated in Git, while figures reference `../figures_en/*.png`.
"""
from __future__ import annotations

from collections import Counter
import hashlib
from html import escape, unescape
import os
from pathlib import Path
import re
import sys
import tempfile
from urllib.parse import unquote, urlsplit, urlunsplit

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import build_site
from scripts.build_markdown_helpers import (
    map_table_cell_text,
    render_markdown,
    validate_url_schemes,
)
from scripts.chapter_identity import LEGACY_HTML_ROUTES, canonical_source, inject_exercise_aliases
from scripts.heading_aliases import has_historical_sequential_aliases, historical_aliases

SRC_EN = ROOT / "chapters_en"
SRC_ZH = ROOT / "chapters"
RESEARCH_EN = ROOT / "codes" / "chapters" / "ch00" / "research_en"
RESEARCH_ZH = ROOT / "codes" / "chapters" / "ch00" / "research"
FIGURES_EN = ROOT / "figures_en"
OUT_EN = ROOT / "site_en"
COURSE_HOME_EN_ROOT = ROOT / "scripts" / "course_home_en"
COURSE_HOME_FILES = ("index.html", "home.css", "home.js")
HOME_FNAME = "00_overview.md"

CHAPTERS_EN = [
    ("00_math-foundations.md", "Chapter 0 · Mathematical & Signal Analysis Foundations"),
    ("01_problem-definition.md", "Chapter 1 · Problem Formulation & Binaural Cues"),
    ("02_basics-signal-model.md", "Chapter 2 · Acoustic Propagation & Array Signal Model"),
    ("03_array-geometry.md", "Chapter 3 · Microphone Array Geometry"),
    ("04_doa-estimation.md", "Chapter 4 · Sound Source Localization (DOA)"),
    ("05_beamforming.md", "Chapter 5 · Beamforming"),
    ("06_aec.md", "Chapter 6 · Acoustic Echo Cancellation (AEC)"),
    ("07_wpe-dereverberation.md", "Chapter 7 · Dereverberation (WPE)"),
    ("08_speech-separation.md", "Chapter 8 · Multichannel Speech Separation"),
    ("09_source-tracking.md", "Chapter 9 · Sound Source Tracking"),
    ("10_engineering-practice.md", "Chapter 10 · Real-Time Systems, Evaluation & Deployment"),
    ("11_selection-guide.md", "Chapter 11 · System Synthesis & Selection Guide"),
    ("12_acoustic-imaging.md", "Chapter 12 · Acoustic Imaging & Noise Diagnostics"),
    ("13_distributed-enhancement.md", "Chapter 13 · Distributed Wireless Microphone Arrays"),
    ("14_spatial-wearable-arrays.md", "Chapter 14 · Spatial Sound Fields & Wearable Arrays"),
    ("15_target-extraction-diarization.md", "Chapter 15 · Conditioned Target Extraction & Diarization"),
    ("16_appendix-symbols-math.md", "Chapter 16 · Appendix A: Symbols, Glossary & Math Lookup"),
    ("17_appendix-guide.md", "Chapter 17 · Appendix B: Literature Map, Triage & Exercises"),
]

RESEARCH_EN_PAGES = [
    ("README.md", "Upstream Source Manuals"),
    ("01_spatial_and_tracking.md", "Spatial Processing & Tracking Sources"),
    ("02_aec_wpe_separation.md", "AEC, WPE & Separation Sources"),
    ("03_industrial_deployment.md", "Industrial Audio Deployment Manual"),
    ("04_source_reproduction.md", "Upstream Acquisition & Reproduction"),
    ("05_exercises_and_audio.md", "Chapter Code Exercises & Audio Suites"),
]

LABEL_BY_FNAME_EN = {fname: label for fname, label in CHAPTERS_EN}
LABEL_BY_FNAME_EN[HOME_FNAME] = "🏠 Course Overview (Home)"

CSS_EN = build_site.CSS.replace(
    'font-family:-apple-system,"PingFang SC","Hiragino Sans GB","Microsoft YaHei",sans-serif',
    'font-family:system-ui,-apple-system,"Segoe UI",Roboto,"Helvetica Neue",Arial,sans-serif',
)

PAGE_EN = """<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title} · Microphone Array Signal Processing</title>
<meta name="source-digest" content="{source_digest}"><style>{css}</style>
<script>
{inline_layout_script}
window.MathJax = {{tex: {{inlineMath: [['$', '$'], ['\\\\(', '\\\\)']], displayMath: [['$$', '$$']]}}, startup: {{pageReady: () => MathJax.startup.defaultPageReady().then(() => document.fonts.ready).then(() => ArrayTutorialLayout.apply())}}}};
</script>
<script defer src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-mml-chtml.js"
 onerror="document.getElementById('offnote').style.display='block';document.getElementById('offnote').textContent='MathJax script failed to load: displaying raw TeX source.';"></script>
</head><body id="top"><a class="skip-link" href="#main-content">Skip to main content</a>
<header class="topbar"><span><a href="{home_href}">🏠 Home</a> &nbsp;/&nbsp; {crumb}</span><a class="lang-switch" href="{zh_url}" data-lang-peer="{zh_peer}" lang="zh-CN" hreflang="zh-CN">中文</a></header>
<div class="wrap"><nav class="side" aria-label="Course contents">{sidebar}</nav>
<main class="main" id="main-content" tabindex="-1"><div class="offline-note" id="offnote">Offline mode: displaying raw TeX formulas; prose and audio remain unaffected.</div>
<details class="toc-mobile"><summary>On this page</summary><nav aria-label="On this page">{toc}</nav></details>
{body}{pn}
<footer class="foot">Microphone Array Signal Processing Tutorial · English static edition built by scripts/build_site_en.py</footer>
</main></div><a class="topbtn" href="#top" title="Back to top" aria-label="Back to top">↑</a>
<script>if(!navigator.onLine)document.getElementById('offnote').style.display='block';for(const a of document.querySelectorAll('a[data-lang-peer]'))if(location.protocol==='file:'||location.hostname==='localhost'||location.hostname==='127.0.0.1')a.href=a.dataset.langPeer;</script>
</body></html>
"""


def source_outputs_en() -> dict[Path, str]:
    outputs = {(SRC_EN / HOME_FNAME).resolve(): "index.html"}
    for fname, _ in CHAPTERS_EN:
        outputs[(SRC_EN / fname).resolve()] = fname.replace(".md", ".html")
    for old, current in LEGACY_HTML_ROUTES.items():
        outputs[(SRC_EN / old.replace(".html", ".md")).resolve()] = old
    for fname, _ in RESEARCH_EN_PAGES:
        html_name = "index.html" if fname == "README.md" else fname.replace(".md", ".html")
        outputs[(RESEARCH_EN / fname).resolve()] = f"research/{html_name}"
    return outputs


def expected_html_set_en() -> set[str]:
    result = {"index.html"}
    result.update(fname.replace(".md", ".html") for fname, _ in CHAPTERS_EN)
    result.update(LEGACY_HTML_ROUTES.keys())
    for fname, _ in RESEARCH_EN_PAGES:
        html_name = "index.html" if fname == "README.md" else fname.replace(".md", ".html")
        result.add(f"research/{html_name}")
    return result


def source_digest_en() -> str:
    from scripts import publication_dependencies_en as dependencies

    paths = dependencies.site_sources(
        ROOT, SRC_EN, RESEARCH_EN, [name for name, _ in RESEARCH_EN_PAGES])
    return dependencies.source_digest(ROOT, paths)


def paired_zh_heading_aliases(source_path: Path) -> list[list[str]]:
    """Return per-heading-index Chinese primary IDs and slugs so all `#sec-u-*` & `#slug` links work."""
    is_research = source_path.resolve().parent == RESEARCH_EN.resolve()
    zh_path = (RESEARCH_ZH if is_research else SRC_ZH) / source_path.name
    if not zh_path.is_file():
        return []
    zh_md = zh_path.read_text(encoding="utf-8")
    zh_records = build_site.heading_records(build_site.parse_headings(zh_md))
    zh_slugs: Counter[str] = Counter()
    aliases_per_idx: list[list[str]] = []
    for _lvl, zh_text, zh_primary, _zh_legacy in zh_records:
        ids = [zh_primary]
        if is_research:
            slug = re.sub(r"\s+", "-", re.sub(r"[^\w\s-]", "", zh_text.lower()))
            zh_slugs[slug] += 1
            if zh_slugs[slug] > 1:
                slug += f"-{zh_slugs[slug] - 1}"
            if slug and slug not in ids:
                ids.append(slug)
        aliases_per_idx.append(ids)
    return aliases_per_idx


def rewrite_images_en(html: str, source_path: Path, current: str) -> str:
    def replace(match):
        parsed, target = build_site.local_link_target(unescape(match.group(3)), source_path)
        if target == (build_site.ROOM_AUDIO_ROOT / "ROOM_RESULTS.png").resolve():
            output = "../site/room_audio/ROOM_RESULTS.png"
        elif (
            target is not None
            and target.parent in {(FIGURES_EN).resolve(), (ROOT / "figures").resolve()}
            and target.suffix.lower() == ".png"
        ):
            output = "../figures_en/" + target.name
        else:
            return match.group(0)
        relative = os.path.relpath(output, Path(current).parent).replace(os.sep, "/")
        value = urlunsplit(("", "", relative, parsed.query, parsed.fragment))
        return match.group(1) + match.group(2) + escape(value, quote=True) + match.group(2)

    return re.sub(r'(<img\b[^>]*?\s+src\s*=\s*)(["\x27])(.*?)\2', replace, html, flags=re.S)


AUDIO_DIR_MAP = [
    (build_site.REAL_AUDIO_ROOT, build_site.REAL_AUDIO_FILES, "real_audio"),
    (build_site.ROOM_AUDIO_ROOT, None, "room_audio"),
    (build_site.TRACKING_AUDIO_ROOT, set(build_site.TRACKING_AUDIO_WAVS) | {"MANIFEST.json"}, "tracking_audio"),
    (build_site.MOVING_AUDIO_ROOT, set(build_site.MOVING_AUDIO_WAVS) | {"MANIFEST.json"}, "moving_audio"),
    (build_site.GSS_AUDIO_ROOT, set(build_site.GSS_AUDIO_WAVS) | {"MANIFEST.json", "STATE.npz"}, "gss_audio"),
    (build_site.BINAURAL_AUDIO_ROOT, build_site.BINAURAL_AUDIO_WAVS | {"MANIFEST.json"}, "binaural_audio"),
    (build_site.SPECTRAL_AUDIO_ROOT, build_site.SPECTRAL_AUDIO_WAVS | {"MANIFEST.json"}, "spectral_audio"),
    (build_site.STFT_AUDIO_ROOT, build_site.STFT_AUDIO_WAVS | {"MANIFEST.json"}, "stft_audio"),
    (build_site.SWEEP_AUDIO_ROOT, build_site.SWEEP_AUDIO_WAVS | {"MANIFEST.json"}, "sweep_audio"),
    (build_site.BASELINE_AUDIO_ROOT, build_site.BASELINE_AUDIO_WAVS | {"MANIFEST.json"}, "baseline_audio"),
    (build_site.REFLECTION_AUDIO_ROOT, build_site.REFLECTION_AUDIO_WAVS | {"MANIFEST.json"}, "reflection_audio"),
    (build_site.BINAURAL_BEAMFORMING_ROOT, build_site.BINAURAL_BEAMFORMING_WAVS | {"MANIFEST.json"}, "binaural_beamforming_audio"),
    (build_site.PHASE_AUDIO_ROOT, build_site.PHASE_AUDIO_WAVS | {"MANIFEST.json"}, "phase_audio"),
    (build_site.GEOMETRY_AUDIO_ROOT, set(build_site.GEOMETRY_AUDIO_WAVS) | {"MANIFEST.json"}, "geometry_audio"),
    (build_site.DERIVATIVE_AUDIO_ROOT, set(build_site.DERIVATIVE_AUDIO_WAVS) | {"MANIFEST.json"}, "derivative_audio"),
    (build_site.APA_AUDIO_ROOT, build_site.APA_AUDIO_WAVS | {"MANIFEST.json"}, "apa_audio"),
    (build_site.REFERENCE_AUDIO_ROOT, build_site.REFERENCE_AUDIO_WAVS | {"MANIFEST.json"}, "reference_audio"),
    (build_site.DELAY_AUDIO_ROOT, build_site.DELAY_AUDIO_WAVS | {"MANIFEST.json"}, "delay_audio"),
    (build_site.CSS_AUDIO_ROOT, build_site.CSS_AUDIO_WAVS | {"MANIFEST.json"}, "css_audio"),
    (build_site.MINT_AUDIO_ROOT, set(build_site.MINT_AUDIO_WAVS) | {"MANIFEST.json"}, "mint_audio"),
    (build_site.SCENARIO_AUDIO_ROOT, build_site.SCENARIO_AUDIO_WAVS | {"MANIFEST.json"}, "scenario_audio"),
    (build_site.RESPONSE_AUDIO_ROOT, build_site.RESPONSE_AUDIO_WAVS | {"MANIFEST.json"}, "response_audio"),
    (build_site.DISTRIBUTED_AUDIO_ROOT, build_site.DISTRIBUTED_AUDIO_WAVS | {"MANIFEST.json"}, "distributed_audio"),
    (build_site.WEARABLE_AUDIO_ROOT, build_site.WEARABLE_AUDIO_WAVS | {"MANIFEST.json"}, "wearable_audio"),
    (build_site.TSE_AUDIO_ROOT, build_site.TSE_AUDIO_WAVS | {"MANIFEST.json"}, "tse_audio"),
    (build_site.IMAGING_AUDIO_ROOT, build_site.IMAGING_AUDIO_WAVS | {"MANIFEST.json"}, "imaging_audio"),
    (build_site.WEIGHTED_AUDIO_ROOT, build_site.WEIGHTED_AUDIO_WAVS | {"MANIFEST.json"}, "weighted_audio"),
    (build_site.NOISE_AUDIO_ROOT, build_site.NOISE_AUDIO_WAVS | {"MANIFEST.json"}, "noise_audio"),
    (build_site.CHANNEL_AUDIO_ROOT, build_site.CHANNEL_AUDIO_WAVS | {"MANIFEST.json"}, "channel_audio"),
    (build_site.MASK_AUDIO_ROOT, build_site.MASK_AUDIO_WAVS | {"MANIFEST.json"}, "mask_audio"),
    (build_site.FOCUS_AUDIO_ROOT, set(build_site.FOCUS_AUDIO_WAVS) | {"MANIFEST.json"}, "focus_audio"),
]


def rewrite_site_links_en(html: str, source_path: Path) -> str:
    outputs = source_outputs_en()
    # Map accidental links targeting Chinese chapters/ or research/ to their English equivalents
    zh_to_en = {(SRC_ZH / fname).resolve(): outputs[(SRC_EN / fname).resolve()]
                for fname in [HOME_FNAME, *(f for f, _ in CHAPTERS_EN)]}
    for old, _cur in LEGACY_HTML_ROUTES.items():
        zh_to_en[(SRC_ZH / old.replace(".html", ".md")).resolve()] = old
    for fname, _ in RESEARCH_EN_PAGES:
        zh_to_en[(RESEARCH_ZH / fname).resolve()] = outputs[(RESEARCH_EN / fname).resolve()]

    current = outputs.get(Path(source_path).resolve(), "index.html")
    main_audio = build_site.main_audio_sources()

    def transform(href: str) -> str:
        parsed, target = build_site.local_link_target(href, source_path)
        if target is None:
            return href
        if target.parent == SRC_EN.resolve():
            canon_target = (SRC_EN / canonical_source(target.name)).resolve()
            if canon_target in outputs:
                relative = os.path.relpath(outputs[canon_target], Path(current).parent).replace(os.sep, "/")
                return urlunsplit(("", "", relative, parsed.query, parsed.fragment))
        if target in outputs:
            relative = os.path.relpath(outputs[target], Path(current).parent).replace(os.sep, "/")
            return urlunsplit(("", "", relative, parsed.query, parsed.fragment))
        if target in zh_to_en:
            relative = os.path.relpath(zh_to_en[target], Path(current).parent).replace(os.sep, "/")
            return urlunsplit(("", "", relative, parsed.query, parsed.fragment))
        if target.parent in {FIGURES_EN.resolve(), (ROOT / "figures").resolve()} and target.suffix.lower() == ".png":
            relative = os.path.relpath("../figures_en/" + target.name, Path(current).parent).replace(os.sep, "/")
            return urlunsplit(("", "", relative, parsed.query, parsed.fragment))
        if target in main_audio:
            relative = os.path.relpath("../site/audio/" + main_audio[target], Path(current).parent).replace(os.sep, "/")
            return urlunsplit(("", "", relative, parsed.query, parsed.fragment))
        for root_dir, allowed_set, folder_name in AUDIO_DIR_MAP:
            if target.parent == root_dir.resolve() and (allowed_set is None or target.name in allowed_set):
                relative = os.path.relpath(f"../site/{folder_name}/{target.name}", Path(current).parent).replace(os.sep, "/")
                return urlunsplit(("", "", relative, parsed.query, parsed.fragment))
        return build_site.repository_url(parsed, target)

    html = build_site.rewrite_href_targets(html, transform)
    html = rewrite_images_en(html, source_path, current)

    audio_folders = "|".join(["audio"] + [folder for _, _, folder in AUDIO_DIR_MAP])
    audio_re = re.compile(rf"^(?:\.\./)+site/(?:{audio_folders})/[a-z0-9_]+\.wav$")

    def player(match):
        href, label = unescape(match.group(1)), match.group(2)
        parsed = urlsplit(href)
        if parsed.path.endswith("real_audio/demand_nriver_16ch_10s.wav"):
            return match.group(0)
        if parsed.scheme or parsed.query or parsed.fragment or not audio_re.fullmatch(parsed.path):
            return match.group(0)
        safe_href = escape(href, quote=True)
        safe_label = escape(re.sub(r"<[^>]+>", "", unescape(label)), quote=True)
        return (
            match.group(0)
            + f'<audio class="audio-sample" controls preload="none" '
            f'aria-label="{safe_label}" src="{safe_href}">Use the WAV link above to download.</audio>'
        )

    return re.sub(r'<a href="([^"]+)">(.*?)</a>', player, html, flags=re.S)


def render_en(md_text: str, source_path: Path) -> tuple[str, int]:
    is_research = source_path.resolve().parent == RESEARCH_EN.resolve()
    records = build_site.heading_records(build_site.parse_headings(md_text))
    zh_aliases = paired_zh_heading_aliases(source_path)
    html = render_markdown(md_text)
    validate_url_schemes(html)
    counter = [0]
    research_slugs: Counter[str] = Counter()
    used_ids: set[str] = {primary for _l, _t, primary, _leg in records}

    def repl(m):
        idx = counter[0]
        counter[0] += 1
        tag, inner = m.group(1), m.group(2)
        _level, text, primary, legacy = records[idx]
        alias_ids: list[str] = []
        if primary != legacy and not has_historical_sequential_aliases(source_path.name):
            alias_ids.append(legacy)
        for old in historical_aliases(source_path.name, primary):
            if old != primary:
                alias_ids.append(old)
        if idx < len(zh_aliases):
            zh_primary = zh_aliases[idx][0]
            for old in historical_aliases(source_path.name, zh_primary):
                if old != primary:
                    alias_ids.append(old)
            for zh_id in zh_aliases[idx]:
                if zh_id != primary:
                    alias_ids.append(zh_id)
        if is_research:
            slug = re.sub(r"\s+", "-", re.sub(r"[^\w\s-]", "", text.lower()))
            research_slugs[slug] += 1
            if research_slugs[slug] > 1:
                slug += f"-{research_slugs[slug] - 1}"
            if slug and slug != primary:
                alias_ids.append(slug)

        alias_spans = []
        for aid in alias_ids:
            if aid and aid not in used_ids:
                used_ids.add(aid)
                alias_spans.append(f'<span id="{escape(aid, quote=True)}" class="anchor-alias" aria-hidden="true"></span>')
        return f'{"".join(alias_spans)}<{tag} id="{primary}">{inner}</{tag}>'

    html = re.sub(r"<(h[1-4])>(.*?)</\1>", repl, html, flags=re.S)
    html = inject_exercise_aliases(html, source_path.name)
    html = re.sub(r'<th(?=[\s>])(?![^>]*\bscope=)([^>]*)>', r'<th scope="col"\1>', html, flags=re.S)
    if source_path.resolve() == (RESEARCH_EN / "02_aec_wpe_separation.md").resolve():
        html = re.sub(
            r'<a id="(bss|neural|tflocoformer|separation-upstream-audit|mffca)"></a>',
            r'<a id="\1" class="separation-topic-anchor"></a>',
            html,
        )
    if source_path.resolve() == (RESEARCH_EN / "04_source_reproduction.md").resolve():
        html = re.sub(
            r'<a id="(distributed-reproduction|historical-source-bindings)"></a>',
            r'<a id="\1" class="source-reproduction-anchor"></a>',
            html,
        )
    html = re.sub(
        r"<table([^>]*)>(.*?)</table>",
        r'<div class="table-scroll" tabindex="0" role="region" aria-label="Scrollable data table"><table\1>\2</table></div>',
        html,
        flags=re.S,
    )
    html = rewrite_site_links_en(html, source_path)

    def nav_text(m):
        href, text = m.group(1), m.group(2).strip()
        base = href.split("#")[0].split("/")[-1]
        if base == "index.html" and text.endswith(".md"):
            return f'<a href="{href}">🏠 Course Overview</a>'
        for fname, label in LABEL_BY_FNAME_EN.items():
            if base == fname.replace(".md", ".html") and text == fname:
                return f'<a href="{href}">{label}</a>'
        return m.group(0)

    html = re.sub(r'<code>(<a href="[^"]+">[^<>]+</a>)</code>', r"\1", html)
    html = re.sub(r'<a href="([^"]+)">(?:<code>)?([^<>]*\.md)(?:</code>)?</a>', nav_text, html)
    html = re.sub(
        r'<img ([^>]*?)src="([^"]+)"',
        r'''<img \1src="\2" onerror="this.alt+=' (figure missing)'"''',
        html,
    )

    def bare_link(m):
        pre, fname = m.group(1), m.group(2)
        prefix = "../" if is_research else ""
        if fname == HOME_FNAME:
            return f'{pre}<a href="{prefix}index.html">{fname}</a>'
        if fname in LABEL_BY_FNAME_EN:
            return f'{pre}<a href="{prefix}{fname.replace(".md", ".html")}">{fname}</a>'
        return m.group(0)

    html = map_table_cell_text(
        html, lambda text: re.sub(r'(^|[\s>(])((?:0\d|1\d)_[^<\s)"]+\.md)', bare_link, text)
    )
    return html, counter[0]


def sidebar_en(current: str | None, heads) -> str:
    parts = []
    if current is None:
        parts.append('<div class="chap cur" aria-current="page">🏠 Course Overview (Home)</div>')
        lst, _ = build_site.sub_list("index.html", heads)
        parts.append(lst)
    else:
        parts.append('<div class="chap"><a href="index.html">🏠 Course Overview (Home)</a></div>')
    for fname, label in CHAPTERS_EN:
        html_name = fname.replace(".md", ".html")
        if fname == current:
            parts.append(f'<div class="chap cur" aria-current="page">{label}</div>')
            lst, _ = build_site.sub_list(html_name, heads, include_level1=True)
            parts.append(lst)
        else:
            parts.append(f'<div class="chap"><a href="{html_name}">{label}</a></div>')
    parts.append('<div class="chap"><a href="research/index.html">Upstream Source Manuals</a></div>')
    return "\n".join(parts)


def research_sidebar_en(current: str, heads) -> str:
    parts = ['<div class="chap"><a href="../index.html">🏠 Course Overview (Home)</a></div>']
    for fname, label in RESEARCH_EN_PAGES:
        html_name = "index.html" if fname == "README.md" else fname.replace(".md", ".html")
        if fname == current:
            parts.append(f'<div class="chap cur" aria-current="page">{label}</div>')
            parts.append(build_site.sub_list(html_name, heads)[0])
        else:
            parts.append(f'<div class="chap"><a href="{html_name}">{label}</a></div>')
    for fname, label in CHAPTERS_EN:
        parts.append(f'<div class="chap"><a href="../{fname.replace(".md", ".html")}">{label}</a></div>')
    return "\n".join(parts)


def with_course_landing_en(page: str) -> str:
    template, css, script = (
        (COURSE_HOME_EN_ROOT / name).read_text(encoding="utf-8") for name in COURSE_HOME_FILES
    )
    if "</header>\n" not in template or not template.endswith("</div>\n"):
        raise ValueError("English course home template missing </header> or closing </div>")
    template_main = template.replace(
        "</header>\n",
        '</header>\n</div>\n<main id="course-page-main">\n<div class="course">\n',
        1,
    )
    return (
        page.replace("</head>", "<style>" + css + "</style></head>", 1)
        .replace('<body id="top">', '<body id="top" class="course-page">', 1)
        .replace('href="#main-content">Skip to main content', 'href="#course-title">Skip to course content', 1)
        .replace(
            '<header class="topbar"><span><a href="index.html">🏠 Home</a> &nbsp;/&nbsp; Course Overview</span><a class="lang-switch" href="https://nanless.github.io/microphone-array-tutorial/site/index.html" data-lang-peer="../site/index.html" lang="zh-CN" hreflang="zh-CN">中文</a></header>\n',
            "",
            1,
        )
        .replace('<div class="wrap">', template_main + '<div class="wrap">', 1)
        .replace('<main class="main" id="main-content" tabindex="-1">', '<div class="main" id="main-content" tabindex="-1">', 1)
        .replace('</main></div><a class="topbtn"', '</div></div></main><a class="topbtn"', 1)
        .replace("</body>", "<script>" + script + "</script></body>", 1)
    )


def main() -> None:
    names = [f for f, _ in CHAPTERS_EN]
    expected = expected_html_set_en()
    pages_zh_base = "https://nanless.github.io/microphone-array-tutorial/site/"
    with tempfile.TemporaryDirectory(prefix=".site-en-build-", dir=ROOT) as tmp:
        temp_out = Path(tmp)
        build_digest = source_digest_en()
        inline_layout_script = build_site.INLINE_LAYOUT_PATH.read_text(encoding="utf-8")
        home_md = (SRC_EN / HOME_FNAME).read_text(encoding="utf-8")
        home_heads = build_site.parse_headings(home_md)
        home_html, home_n = render_en(home_md, SRC_EN / HOME_FNAME)
        assert home_n == len(home_heads), f"Home anchors {home_n} vs headings {len(home_heads)}"
        toc, _ = build_site.sub_list("index.html", home_heads)
        home_html = re.sub(r"<h1(\b[^>]*>.*?)</h1>", r"<h2\1</h2>", home_html, count=1, flags=re.S)
        home_page = PAGE_EN.format(
            title="From Fundamentals to Industrial Practice",
            css=CSS_EN,
            crumb="Course Overview",
            home_href="index.html",
            zh_url=pages_zh_base + "index.html",
            zh_peer="../site/index.html",
            source_digest=build_digest,
            inline_layout_script=inline_layout_script,
            sidebar=sidebar_en(None, home_heads),
            toc=toc,
            body=home_html,
            pn="",
        )
        (temp_out / "index.html").write_text(with_course_landing_en(home_page), encoding="utf-8")

        for i, (fname, label) in enumerate(CHAPTERS_EN):
            md = (SRC_EN / fname).read_text(encoding="utf-8")
            heads = build_site.parse_headings(md)
            html_name = fname.replace(".md", ".html")
            body, n = render_en(md, SRC_EN / fname)
            assert n == len(heads), f"{fname}: anchors {n} vs headings {len(heads)}"
            body, heads = build_site.promote_content_headings(body, heads)
            toc, _ = build_site.sub_list(html_name, heads, include_level1=True)
            prev = (
                f'<a href="{names[i - 1].replace(".md", ".html")}">← Previous Chapter</a>'
                if i > 0
                else '<a href="index.html">← Course Overview</a>'
            )
            nxt = (
                f'<a href="{names[i + 1].replace(".md", ".html")}">Next Chapter →</a>'
                if i < len(names) - 1
                else '<a href="index.html">Back to Home →</a>'
            )
            pn = f'<div class="pn"><span>{prev}</span><span>{nxt}</span></div>'
            (temp_out / html_name).write_text(
                PAGE_EN.format(
                    title=label,
                    css=CSS_EN,
                    crumb=label,
                    home_href="index.html",
                    zh_url=pages_zh_base + html_name,
                    zh_peer="../site/" + html_name,
                    source_digest=build_digest,
                    inline_layout_script=inline_layout_script,
                    sidebar=sidebar_en(fname, heads),
                    toc=toc,
                    body=body,
                    pn=pn,
                ),
                encoding="utf-8",
            )

        for old, current in LEGACY_HTML_ROUTES.items():
            legacy = (temp_out / current).read_text(encoding="utf-8")
            legacy = legacy.replace("<head>", f'<head>\n<link rel="canonical" href="{current}">', 1)
            (temp_out / old).write_text(legacy, encoding="utf-8")

        (temp_out / "research").mkdir()
        for fname, label in RESEARCH_EN_PAGES:
            source = RESEARCH_EN / fname
            md = source.read_text(encoding="utf-8")
            heads = build_site.parse_headings(md)
            body, n = render_en(md, source)
            assert n == len(heads), f"research_en/{fname}: anchors {n} vs headings {len(heads)}"
            html_name = "index.html" if fname == "README.md" else fname.replace(".md", ".html")
            (temp_out / "research" / html_name).write_text(
                PAGE_EN.format(
                    title=label,
                    css=CSS_EN,
                    crumb=label,
                    home_href="../index.html",
                    zh_url=pages_zh_base + "research/" + html_name,
                    zh_peer="../../site/research/" + html_name,
                    source_digest=build_digest,
                    sidebar=research_sidebar_en(fname, heads),
                    inline_layout_script=inline_layout_script,
                    toc=build_site.sub_list(html_name, heads)[0],
                    body=body,
                    pn="",
                ),
                encoding="utf-8",
            )

        built = {path.relative_to(temp_out).as_posix() for path in temp_out.rglob("*.html")}
        if built != expected:
            raise SystemExit(f"English site page set mismatch: expected {sorted(expected)}, got {sorted(built)}")

        OUT_EN.mkdir(exist_ok=True)
        (OUT_EN / "research").mkdir(exist_ok=True)
        stale = [path for path in OUT_EN.glob("*.html") if path.name not in expected]
        stale += [
            path
            for path in (OUT_EN / "research").glob("*.html")
            if "research/" + path.name not in expected
        ]
        build_site.publish_files(
            [(temp_out / name, OUT_EN / name) for name in sorted(expected)],
            stale,
            boundary=OUT_EN,
        )
    print("DONE EN", len(expected), "pages")


if __name__ == "__main__":
    main()
