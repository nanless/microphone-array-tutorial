"""Explicit English publication inputs, shared by producers and freshness gates."""
from __future__ import annotations

import hashlib
from pathlib import Path

from scripts.figure_report_contracts import FIGURE_REPORTS


SHARED_SCRIPTS = (
    "build_site.py", "build_markdown_helpers.py", "inline_layout.js",
    "chapter_identity.py", "chapter_numbering.json", "heading_aliases.py",
    "legacy_sequential_anchors.json", "publication_dependencies_en.py",
    "make_figures_en.py", "figure_en_translations.json",
    "make_figures.py", "make_aec_figures.py", "make_beamforming_figures.py",
    "make_reference_figures.py", "make_delay_figures.py", "make_css_figures.py",
    "make_tracking_figures.py", "make_channel_figures.py", "make_selection_figures.py",
    "make_wearable_figures.py", "make_tse_figures.py", "make_math_figures.py",
    "figure_report_contracts.py", "portable_figure_report_checks.py",
    "pdf_input_identity.py",
    "quality_check.py",
)


def common_sources(root: Path, source: Path, research: Path, research_names) -> list[Path]:
    """Include paired Chinese headings used to emit compatible English anchors."""
    chapters = sorted(source.glob("*.md"))
    manuals = [research / name for name in research_names]
    paths = chapters + manuals
    paths += [root / "chapters" / p.name for p in chapters]
    paths += [root / "codes/chapters/ch00/research" / p.name for p in manuals]
    paths += [root / "scripts" / name for name in SHARED_SCRIPTS]
    paths += sorted((root / "figures_en").glob("fig*.png"))
    paths += [root / entry.report_relative for entry in FIGURE_REPORTS]
    paths.append(root / "requirements.txt")
    return paths


def site_sources(root: Path, source: Path, research: Path, research_names) -> list[Path]:
    paths = common_sources(root, source, research, research_names)
    paths.append(root / "scripts/build_site_en.py")
    paths += [root / "scripts/course_home_en" / name
              for name in ("index.html", "home.css", "home.js")]
    paths += [root / "assets/readme/cases/locata" / name for name in
              ("MANIFEST.json", "NOTICE.txt", "before.wav", "after.wav",
               "mvdr.wav", "focus.wav", "levels.svg")]
    paths += [root / "assets/readme/lab" / name
              for name in ("array-lab.html", "single.wav", "left.wav", "right.wav")]
    # English players link to these same chapter assets through the Chinese site.
    for manifest in sorted((root / "codes/chapters").glob("*/*/MANIFEST.json")):
        paths += sorted(p for p in manifest.parent.iterdir() if p.is_file())
    paths += sorted((root / "codes/chapters").glob("*/audio/*.wav"))
    return paths


def pdf_sources(root: Path, source: Path, research: Path, research_names) -> list[Path]:
    paths = common_sources(root, source, research, research_names)
    paths.append(root / "codes/chapters/appendix_b/room_audio/ROOM_RESULTS.png")
    paths += [root / "scripts" / name for name in
              ("build_pdf_en.py", "build_pdf.py", "build_site_en.py", "pdf_math.js")]
    # The website uses a CDN; only PDF rendering consumes these local resources.
    paths += sorted(p for p in (root / "scripts/vendor/mathjax-3.2.2").rglob("*")
                    if p.is_file())
    return paths


def source_digest(root: Path, paths) -> str:
    digest = hashlib.sha256()
    for path in sorted(set(paths)):
        digest.update(path.relative_to(root).as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()[:12]
