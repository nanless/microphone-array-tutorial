# -*- coding: utf-8 -*-
"""Generate all 87 English scientific figures into ``figures_en/`` non-destructively.

Usage (from repository root):
    .venv/bin/python scripts/make_figures_en.py

Design guarantees:
1. Never modifies any of the 12 Chinese figure scripts, ``figures/*.png``, or
   ``codes/chapters/*/reports/*.json`` (all report writes during execution are
   redirected to an isolated temporary directory).
2. Configures clean academic sans-serif English typography (``Helvetica Neue``,
   ``Arial``, ``DejaVu Sans``), with automatic per-artist fallback to
   ``DejaVu Sans`` whenever a string contains mathematical/arrow symbols outside
   ``Helvetica Neue``'s glyph table.
3. Translates every ``matplotlib.text.Text`` artist, axis ``Formatter``
   (including ``FuncFormatter``, ``StrCategoryFormatter``, ``FixedFormatter``),
   and ``Table`` cell via ``scripts/figure_en_translations.json`` (exact lookup
   + regex rules), applies geometry-aware font-size fitting so translated
   strings never overflow flowchart boxes, table cells, subplot titles, or
   figure bounds, and asserts ZERO remaining CJK/fullwidth characters
   (``[\\u4e00-\\u9fff\\u3000-\\u303f\\uff00-\\uffef]``) before and after ``savefig``.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import runpy
import sys
import tempfile

REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
from matplotlib.ft2font import FT2Font
from matplotlib.patches import FancyBboxPatch, Rectangle
import matplotlib.pyplot as plt
import matplotlib.text as mtext
import matplotlib.ticker as mticker

import codes.chapters.ch00.io_contracts as io_contracts

ENGLISH_SANS_FONTS = ["Helvetica Neue", "Arial", "DejaVu Sans"]
CJK_RE = re.compile(r"[\u4e00-\u9fff\u3000-\u303f\uff00-\uffef]")
MATH_SPAN_RE = re.compile(r"\$[^$]*\$")

FIGURE_SCRIPTS = (
    "make_figures.py",
    "make_aec_figures.py",
    "make_beamforming_figures.py",
    "make_channel_figures.py",
    "make_css_figures.py",
    "make_delay_figures.py",
    "make_math_figures.py",
    "make_reference_figures.py",
    "make_selection_figures.py",
    "make_tracking_figures.py",
    "make_tse_figures.py",
    "make_wearable_figures.py",
)

TRANSLATIONS_PATH = Path(__file__).resolve().parent / "figure_en_translations.json"
FIGURES_ZH_DIR = (REPOSITORY_ROOT / "figures").resolve()
FIGURES_EN_DIR = (REPOSITORY_ROOT / "figures_en").resolve()
CODES_DIR = (REPOSITORY_ROOT / "codes").resolve()


def _load_primary_charmap() -> set[int] | None:
    try:
        font_path = fm.findfont(
            fm.FontProperties(family=ENGLISH_SANS_FONTS[0]),
            fallback_to_default=False,
        )
        return set(FT2Font(font_path).get_charmap().keys())
    except Exception:
        return None


PRIMARY_CHARMAP = _load_primary_charmap()


def _needs_unicode_fallback_font(text: str) -> bool:
    if PRIMARY_CHARMAP is None or not text:
        return False
    outside_math = MATH_SPAN_RE.sub("", text)
    for ch in outside_math:
        if ch in ("\n", "\r", "\t"):
            continue
        if ord(ch) > 127 and ord(ch) not in PRIMARY_CHARMAP:
            return True
    return False


class FigureTranslator:
    """Exact + regex translator backed by ``figure_en_translations.json``."""

    def __init__(self, json_path: Path) -> None:
        payload = json.loads(json_path.read_text(encoding="utf-8"))
        self.exact: dict[str, str] = dict(payload["exact"])
        self.regex_rules: list[tuple[re.Pattern[str], str]] = [
            (re.compile(rule["pattern"]), rule["replacement"])
            for rule in payload.get("regex_rules", [])
        ]

    def translate_one(self, text: str) -> str:
        if not text or not CJK_RE.search(text):
            return text
        if text in self.exact:
            return self.exact[text]
        stripped = text.strip()
        if stripped in self.exact:
            lead = text[: len(text) - len(text.lstrip())]
            tail = text[len(text.rstrip()) :]
            return f"{lead}{self.exact[stripped]}{tail}"
        out = text
        for pattern, repl in self.regex_rules:
            out = pattern.sub(repl, out)
        if not CJK_RE.search(out):
            return out
        if "\n" in text:
            lines = [self.translate_one(line) for line in text.split("\n")]
            joined = "\n".join(lines)
            if not CJK_RE.search(joined):
                return joined
        return out


class _TranslatedFormatter(mticker.Formatter):
    """Wraps any Matplotlib Formatter so regenerated tick labels stay English."""

    def __init__(self, inner: mticker.Formatter, translator: FigureTranslator) -> None:
        self._inner = inner
        self._translator = translator
        self.axis = inner.axis

    def set_axis(self, axis) -> None:
        super().set_axis(axis)
        self._inner.set_axis(axis)

    def set_locs(self, locs) -> None:
        self.locs = locs
        self._inner.set_locs(locs)

    def __call__(self, x, pos=None) -> str:
        s = self._inner(x, pos)
        return self._translator.translate_one(s) if isinstance(s, str) else s

    def format_ticks(self, values):
        return [
            self._translator.translate_one(s) if isinstance(s, str) else s
            for s in self._inner.format_ticks(values)
        ]

    def format_data(self, value) -> str:
        s = self._inner.format_data(value)
        return self._translator.translate_one(s) if isinstance(s, str) else s

    def format_data_short(self, value) -> str:
        s = self._inner.format_data_short(value)
        return self._translator.translate_one(s) if isinstance(s, str) else s

    def get_offset(self) -> str:
        s = self._inner.get_offset()
        return self._translator.translate_one(s) if isinstance(s, str) else s


def _set_artist_font(obj: mtext.Text, text: str) -> None:
    families = obj.get_fontfamily()
    if families and any("mono" in f.lower() for f in families):
        return
    if _needs_unicode_fallback_font(text):
        obj.set_fontfamily(["DejaVu Sans"])
    else:
        obj.set_fontfamily(list(ENGLISH_SANS_FONTS))


def _assert_zero_cjk(fig: matplotlib.figure.Figure, fig_name: str, stage: str) -> None:
    offenders: list[str] = []
    for obj in fig.findobj(mtext.Text, include_self=True):
        s = obj.get_text()
        if s and CJK_RE.search(s):
            offenders.append(s)
    for ax in fig.axes:
        for table in ax.tables:
            for cell in table.get_celld().values():
                s = cell.get_text().get_text()
                if s and CJK_RE.search(s):
                    offenders.append(s)
    if offenders:
        unique = list(dict.fromkeys(offenders))
        raise AssertionError(
            f"[{stage}] {fig_name} still contains CJK/fullwidth characters: {unique!r}"
        )


def _autofit_translated_layout(
    fig: matplotlib.figure.Figure,
    was_translated: set[int],
) -> None:
    """Adjust font sizes gently where English strings would overflow containers."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    fig_bbox = fig.bbox

    for ax in fig.axes:
        ax_bbox = ax.get_window_extent(renderer)
        candidate_patches = []
        for patch in ax.patches:
            if isinstance(patch, (Rectangle, FancyBboxPatch)) and patch.get_visible():
                pb = patch.get_window_extent(renderer)
                if pb.width < 18 or pb.height < 12:
                    continue
                if ax.axison and pb.width > 0.90 * ax_bbox.width and pb.height > 0.90 * ax_bbox.height:
                    continue
                if not ax.axison and pb.width > 0.96 * ax_bbox.width and pb.height > 0.94 * ax_bbox.height:
                    continue
                candidate_patches.append(pb)

        for txt in ax.texts:
            s = txt.get_text().strip()
            if not s or not txt.get_visible():
                continue
            tb = txt.get_window_extent(renderer)
            if tb.width <= 0 or tb.height <= 0:
                continue
            tcx = 0.5 * (tb.x0 + tb.x1)
            tcy = 0.5 * (tb.y0 + tb.y1)

            matched_box = None
            best_area = float("inf")
            for pb in candidate_patches:
                if pb.x0 <= tcx <= pb.x1 and pb.y0 <= tcy <= pb.y1:
                    pcx = 0.5 * (pb.x0 + pb.x1)
                    pcy = 0.5 * (pb.y0 + pb.y1)
                    if (
                        abs(tcx - pcx) <= 0.26 * pb.width
                        and abs(tcy - pcy) <= 0.32 * pb.height
                        and (pb.width * pb.height) < best_area
                    ):
                        matched_box = pb
                        best_area = pb.width * pb.height

            if matched_box is not None:
                max_w = 0.90 * matched_box.width
                max_h = 0.86 * matched_box.height
                if tb.width > max_w or tb.height > max_h:
                    scale = min(max_w / tb.width, max_h / tb.height)
                    txt.set_fontsize(max(7.2, txt.get_fontsize() * scale))
            elif not ax.axison and id(txt) in was_translated:
                for pb in candidate_patches:
                    if pb.x0 + 4 <= tb.x0 <= pb.x1 - 24 and pb.y0 <= tcy <= pb.y1 and tb.x1 > pb.x1 - 10:
                        avail_w = max(40.0, pb.x1 - 12.0 - tb.x0)
                        scale = min(1.0, avail_w / tb.width)
                        txt.set_fontsize(max(7.8, txt.get_fontsize() * scale))
                        break
            elif ax.axison and id(txt) in was_translated and ax_bbox.width > 40:
                if tb.x0 >= ax_bbox.x0 - 4 and tb.x1 > ax_bbox.x1 + 6:
                    avail = max(40.0, ax_bbox.x1 - tb.x0 - 4.0)
                    scale = max(0.76, avail / tb.width)
                    txt.set_fontsize(max(8.0, txt.get_fontsize() * scale))
                elif tb.x1 <= ax_bbox.x1 + 4 and tb.x0 < ax_bbox.x0 - 6:
                    avail = max(40.0, tb.x1 - ax_bbox.x0 - 4.0)
                    scale = max(0.76, avail / tb.width)
                    txt.set_fontsize(max(8.0, txt.get_fontsize() * scale))

        for table in ax.tables:
            for cell in table.get_celld().values():
                cb = cell.get_window_extent(renderer)
                ct = cell.get_text()
                if not ct.get_text().strip() or cb.width < 12:
                    continue
                tb = ct.get_window_extent(renderer)
                if tb.width > 0.90 * cb.width:
                    scale = (0.90 * cb.width) / tb.width
                    ct.set_fontsize(max(7.0, ct.get_fontsize() * scale))

        if ax_bbox.width > 35:
            max_title_w = (1.06 if len(fig.axes) > 1 else 1.22) * ax_bbox.width
            for title_obj in (ax.title, ax._left_title, ax._right_title):
                if title_obj.get_text().strip():
                    tb = title_obj.get_window_extent(renderer)
                    if tb.width > max_title_w:
                        scale = max_title_w / tb.width
                        title_obj.set_fontsize(max(8.5, title_obj.get_fontsize() * scale))

        leg = ax.get_legend()
        if leg is not None and ax_bbox.width > 45:
            lb = leg.get_window_extent(renderer)
            if lb.width > 0.80 * ax_bbox.width:
                scale = max(0.76, (0.78 * ax_bbox.width) / lb.width)
                for lt in leg.get_texts():
                    lt.set_fontsize(max(8.0, lt.get_fontsize() * scale))

        if ax.axison:
            xt_items = []
            for t in ax.get_xticklabels():
                if t.get_visible() and t.get_text().strip():
                    tb = t.get_window_extent(renderer)
                    if tb.width > 0:
                        xt_items.append((0.5 * (tb.x0 + tb.x1), 0.5 * (tb.y0 + tb.y1), tb, t))
            xt_items.sort(key=lambda item: item[0])
            worst_ratio = 1.0
            for idx in range(len(xt_items) - 1):
                cx0, cy0, tb0, _ = xt_items[idx]
                cx1, cy1, tb1, _ = xt_items[idx + 1]
                if abs(cy0 - cy1) < 14.0 and tb0.x1 + 5.0 > tb1.x0 and cx1 > cx0 + 4.0:
                    needed = 0.5 * (tb0.width + tb1.width) + 6.0
                    avail = cx1 - cx0
                    worst_ratio = min(worst_ratio, avail / needed)
            if worst_ratio < 1.0:
                cur_fs = xt_items[0][3].get_fontsize()
                new_fs = max(7.2, cur_fs * worst_ratio)
                ax.tick_params(axis="x", labelsize=new_fs)
                for _, _, _, t in xt_items:
                    t.set_fontsize(new_fs)

    if fig_bbox.width > 60:
        max_fig_text_w = 0.95 * fig_bbox.width
        for ftxt in fig.texts:
            if ftxt.get_text().strip():
                tb = ftxt.get_window_extent(renderer)
                if tb.width > max_fig_text_w:
                    scale = max_fig_text_w / tb.width
                    ftxt.set_fontsize(max(8.5, ftxt.get_fontsize() * scale))


def generate_all_english_figures(only_scripts: tuple[str, ...] = FIGURE_SCRIPTS) -> list[Path]:
    translator = FigureTranslator(TRANSLATIONS_PATH)
    FIGURES_EN_DIR.mkdir(parents=True, exist_ok=True)
    saved_files: list[Path] = []

    with tempfile.TemporaryDirectory(prefix="masp_fig_en_reports_") as tmp_str:
        tmpdir = Path(tmp_str).resolve()

        orig_write_text = Path.write_text
        orig_write_bytes = Path.write_bytes
        orig_read_bytes = Path.read_bytes
        orig_chmod = os.chmod
        orig_write_json_report = io_contracts.write_json_report
        orig_savefig = matplotlib.figure.Figure.savefig

        def _is_protected_repo_write(path: Path) -> Path | None:
            rp = path.resolve()
            for protected_root in (CODES_DIR, FIGURES_ZH_DIR):
                try:
                    rel = rp.relative_to(protected_root.parent)
                    redirected = tmpdir / rel
                    redirected.parent.mkdir(parents=True, exist_ok=True)
                    return redirected
                except ValueError:
                    continue
            return None

        def patched_write_text(self: Path, data: str, *args, **kwargs):
            redir = _is_protected_repo_write(self)
            target = redir if redir is not None else self
            return orig_write_text(target, data, *args, **kwargs)

        def patched_write_bytes(self: Path, data: bytes, *args, **kwargs):
            redir = _is_protected_repo_write(self)
            target = redir if redir is not None else self
            return orig_write_bytes(target, data, *args, **kwargs)

        def patched_read_bytes(self: Path) -> bytes:
            if self.exists():
                return orig_read_bytes(self)
            redir = _is_protected_repo_write(self)
            if redir is not None and redir.exists():
                return orig_read_bytes(redir)
            return orig_read_bytes(self)

        def patched_write_json_report(destination, payload, *, forbidden_roots=()):
            dest_path = Path(destination)
            redir = _is_protected_repo_write(dest_path)
            target = redir if redir is not None else dest_path
            target.parent.mkdir(parents=True, exist_ok=True)
            encoded = json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
            orig_write_text(target, encoded, encoding="utf-8")
            return target

        def patched_chmod(path, mode, *args, **kwargs):
            rp = Path(path).resolve()
            if rp.parent == FIGURES_ZH_DIR:
                en_path = FIGURES_EN_DIR / rp.name
                if en_path.exists():
                    return orig_chmod(en_path, mode, *args, **kwargs)
                return None
            try:
                rp.relative_to(CODES_DIR)
                return None
            except ValueError:
                return orig_chmod(path, mode, *args, **kwargs)

        def patched_savefig(self: matplotlib.figure.Figure, fname, *args, **kwargs):
            in_path = Path(fname).resolve()
            if in_path.parent == FIGURES_ZH_DIR:
                out_path = FIGURES_EN_DIR / in_path.name
            elif in_path.parent == FIGURES_EN_DIR:
                out_path = in_path
            else:
                out_path = FIGURES_EN_DIR / in_path.name

            # 1. Materialize lazy ticks while original fonts are active
            self.canvas.draw()

            # 2. Wrap all axis formatters & translate table cells
            for ax in self.axes:
                axes_list = [ax.xaxis, ax.yaxis]
                if hasattr(ax, "zaxis"):
                    axes_list.append(ax.zaxis)
                for axis in axes_list:
                    maj = axis.get_major_formatter()
                    if not isinstance(maj, _TranslatedFormatter):
                        axis.set_major_formatter(_TranslatedFormatter(maj, translator))
                    minf = axis.get_minor_formatter()
                    if not isinstance(minf, _TranslatedFormatter):
                        axis.set_minor_formatter(_TranslatedFormatter(minf, translator))
                for table in ax.tables:
                    for cell in table.get_celld().values():
                        ct = cell.get_text()
                        orig_s = ct.get_text()
                        new_s = translator.translate_one(orig_s)
                        if new_s != orig_s:
                            ct.set_text(new_s)
                        _set_artist_font(ct, new_s)

            # 3. Switch to clean English sans-serif stack & translate every Text artist
            prev_family = list(plt.rcParams["font.family"])
            prev_sans = list(plt.rcParams["font.sans-serif"])
            prev_minus = plt.rcParams["axes.unicode_minus"]
            try:
                plt.rcParams["font.family"] = "sans-serif"
                plt.rcParams["font.sans-serif"] = list(ENGLISH_SANS_FONTS)
                plt.rcParams["axes.unicode_minus"] = True

                self.canvas.draw()

                was_translated: set[int] = set()
                for obj in self.findobj(mtext.Text, include_self=True):
                    orig_s = obj.get_text()
                    new_s = translator.translate_one(orig_s)
                    if new_s != orig_s:
                        obj.set_text(new_s)
                        was_translated.add(id(obj))
                    _set_artist_font(obj, new_s)

                _assert_zero_cjk(self, out_path.name, "pre-save")
                _autofit_translated_layout(self, was_translated)

                # Re-sync tick label fonts after draw regenerated any tick labels
                for obj in self.findobj(mtext.Text, include_self=True):
                    _set_artist_font(obj, obj.get_text())

                out_path.parent.mkdir(parents=True, exist_ok=True)
                result = orig_savefig(self, out_path, *args, **kwargs)
                _assert_zero_cjk(self, out_path.name, "post-save")
                saved_files.append(out_path)
                return result
            finally:
                plt.rcParams["font.family"] = prev_family
                plt.rcParams["font.sans-serif"] = prev_sans
                plt.rcParams["axes.unicode_minus"] = prev_minus

        Path.write_text = patched_write_text
        Path.write_bytes = patched_write_bytes
        Path.read_bytes = patched_read_bytes
        os.chmod = patched_chmod
        io_contracts.write_json_report = patched_write_json_report
        matplotlib.figure.Figure.savefig = patched_savefig

        orig_argv = sys.argv[:]
        try:
            for script_name in only_scripts:
                script_path = REPOSITORY_ROOT / "scripts" / script_name
                sys.argv = [str(script_path)]
                print(f"[make_figures_en] Executing {script_name} ...", flush=True)
                try:
                    runpy.run_path(str(script_path), run_name="__main__")
                except SystemExit as exc:
                    if exc.code not in (0, None):
                        raise
        finally:
            sys.argv = orig_argv
            Path.write_text = orig_write_text
            Path.write_bytes = orig_write_bytes
            Path.read_bytes = orig_read_bytes
            os.chmod = orig_chmod
            io_contracts.write_json_report = orig_write_json_report
            matplotlib.figure.Figure.savefig = orig_savefig

    return saved_files


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scripts",
        nargs="*",
        default=list(FIGURE_SCRIPTS),
        help="Subset of figure scripts to run (defaults to all 12 scripts = 87 figures).",
    )
    args = parser.parse_args()
    saved = generate_all_english_figures(tuple(args.scripts))
    unique_names = sorted({p.name for p in saved})
    print(f"[make_figures_en] Completed: {len(unique_names)} unique figures written to {FIGURES_EN_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
