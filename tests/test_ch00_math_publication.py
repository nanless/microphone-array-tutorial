"""Publication, figure-provenance, and markdown-contract tests for Chapter 0."""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import re
import shutil
import stat
import tempfile
import unittest

from PIL import Image

from codes.chapters.ch00 import chapter00_exercises as ex00
from codes.chapters.ch00.io_contracts import strict_json_loads, validate_parent_chain
from scripts import quality_check as qc

ROOT = Path(__file__).resolve().parents[1]

CH00_FIG_SPECS = (
    (84, "fig84_complex_phasor_fourier.png", "figure84_complex_phasor_fourier.json"),
    (85, "fig85_subspace_projection_svd.png", "figure85_subspace_projection_svd.json"),
    (86, "fig86_hermitian_rayleigh_gev.png", "figure86_hermitian_rayleigh_gev.json"),
    (87, "fig87_wirtinger_probability_bayesian.png", "figure87_wirtinger_probability_bayesian.json"),
)

EXPECTED_CH00_PLOT_SOURCES = {
    "codes/chapters/ch00/core/chapter00_math.py",
    "codes/chapters/ch00/io_contracts.py",
    "codes/chapters/ch02/core/conventions.py",
    "scripts/make_figures.py",
    "scripts/make_math_figures.py",
}

FORBIDDEN_CALQUES = (
    "单频别名",
    "空间别名",
    "语音缺席",
    "先验缺席",
    "静态合同",
    "数学控制",
    "代码控制",
    "手算夹具",
    "有限支持",
)


def verify_ch00_figures_and_reports(root: Path) -> list[str]:
    """Validate fig84..87 PNGs and figure84..87 JSON reports under root."""
    errors: list[str] = []
    try:
        for num, png_name, json_name in CH00_FIG_SPECS:
            report_path = validate_parent_chain(root / "codes/chapters/ch00/reports" / json_name)
            png_path = validate_parent_chain(root / "figures" / png_name)
            report_bytes = report_path.read_bytes()
            report = strict_json_loads(report_bytes)

            if (
                type(report.get("figure")) is not int
                or report["figure"] != num
                or not isinstance(report.get("exercise_ids"), list)
                or not report["exercise_ids"]
                or any(type(x) is not str or x not in ex00.EXERCISE_IDS for x in report["exercise_ids"])
            ):
                errors.append(f"figure{num} header or exercise_ids contract invalid")
                continue

            sources = report.get("source_sha256")
            if not isinstance(sources, dict) or set(sources) != EXPECTED_CH00_PLOT_SOURCES:
                errors.append(f"figure{num} source_sha256 keyset mismatch")
                continue

            for rel, digest in sources.items():
                actual_sha = hashlib.sha256(validate_parent_chain(root / rel).read_bytes()).hexdigest()
                if digest != actual_sha:
                    errors.append(f"figure{num} source SHA mismatch for {rel}")

            anchors = report.get("anchors")
            if not isinstance(anchors, dict) or not anchors:
                errors.append(f"figure{num} anchors missing")
                continue
            for key, val in anchors.items():
                if isinstance(val, list):
                    if any(type(item) not in (int, float) for item in val):
                        errors.append(f"figure{num} anchor {key} has invalid element type")
                elif type(val) not in (int, float):
                    errors.append(f"figure{num} anchor {key} must be int/float, not {type(val).__name__}")

            with Image.open(png_path) as img:
                width, height = img.size
                metadata = dict(img.info)
            if width < 800 or height < 300:
                errors.append(f"figure{num} resolution too small: {width}x{height}")

            expected_report_sha = hashlib.sha256(report_bytes).hexdigest()
            if (
                metadata.get("SourceScript") != "scripts/make_math_figures.py"
                or metadata.get("SourceScriptDigest") != sources["scripts/make_math_figures.py"]
                or metadata.get("NumericalReportDigest") != expected_report_sha
                or strict_json_loads(metadata.get("GeneratorInputs", "{}").encode("utf-8")) != sources
            ):
                errors.append(f"figure{num} PNG provenance metadata mismatch")
    except (OSError, ValueError, KeyError, TypeError) as exc:
        errors.append(f"Chapter 0 figure verification failed: {exc}")
    return errors


class Chapter00MarkdownAndPublicationTests(unittest.TestCase):
    def setUp(self):
        self.md_path = ROOT / "chapters/00_math-foundations.md"
        self.text = self.md_path.read_text(encoding="utf-8")

    def test_markdown_structure_banners_headings_and_line_budget(self):
        lines = self.text.splitlines()
        self.assertGreaterEqual(len(lines), 1150)
        self.assertLessEqual(len(lines), 1750)
        self.assertTrue(any(line.startswith("> 本篇是第 0 章（正文与扩展专题共 16 章，即第 0～15 章") for line in lines[:5]))
        self.assertIn("导读首页：[00_overview.md](00_overview.md) ｜ 下一篇：[01_problem-definition.md](01_problem-definition.md)", self.text)
        self.assertEqual(lines[-1].strip(), "> 📄 本篇信息：配图 4 张 ｜ [回首页](00_overview.md)")

        headings = qc.markdown_headings(self.text)
        h2 = [t for lvl, t in headings if lvl == 2]
        h3 = [t for lvl, t in headings if lvl == 3]
        self.assertEqual(h2, ["0. 全书数学基础与信号分析预备"])
        self.assertEqual(len(h3), 9)
        for idx in range(1, 10):
            self.assertTrue(h3[idx - 1].startswith(f"0.{idx} "), msg=h3[idx - 1])

        found_ex = qc.exercise_ids_from_markdown(self.text)
        self.assertEqual(found_ex, list(ex00.EXERCISE_IDS))
        for idx in range(1, 17):
            self.assertIn(f'<a id="e00-{idx:02d}"></a>', self.text)

    def test_formulas_cross_links_tables_and_style_guards(self):
        self.assertEqual(qc.source_control_character_issues(self.text), [])
        prose = qc.strip_fenced_code(self.text)
        style_text = re.sub(r"`[^`]*`", "", prose)
        style_text = re.sub(r"https?://[^\s)]+", "", style_text)
        self.assertIsNone(qc.NUMERIC_ASCII_RANGE.search(style_text))
        self.assertIsNone(qc.EDITING_MARKERS.search(self.text))
        self.assertIsNone(re.search(r"\^[T]|\\mathsf\{T\}", prose))
        for bad in FORBIDDEN_CALQUES:
            self.assertNotIn(bad, self.text, msg=f"forbidden calque found: {bad}")

        # Lookup table header & zero raw '|' inside inline $...$ in table rows
        expected_table_header = "| 章节 | 核心声学与信号处理算法 | 依赖的本章数学预备节与核心公式 | 常见数学误区与数值陷阱 |"
        self.assertIn(expected_table_header, self.text)
        for line_no, line in enumerate(prose.splitlines(), 1):
            if line.lstrip().startswith("|"):
                for inline_math in re.findall(r"\$([^$]+)\$", line):
                    self.assertNotIn(
                        "|",
                        inline_math,
                        msg=f"unescaped '|' inside table inline math at line {line_no}: {inline_math}",
                    )

        # Sequential \tag{0-1}..\tag{0-N} and cross-document formula/section references
        docs = {p.name: p.read_text(encoding="utf-8") for p in sorted((ROOT / "chapters").glob("*.md"))}
        ch00_formula_issues = [
            msg for msg in qc.formula_semantic_issues(docs) if "00_math-foundations.md" in msg
        ]
        self.assertEqual(ch00_formula_issues, [])
        ch00_sec_issues = [
            msg for msg in qc.section_reference_issues(docs) if "00_math-foundations.md" in msg
        ]
        self.assertEqual(ch00_sec_issues, [])

        fig_refs = qc.extract_figure_references(self.text)
        self.assertEqual([num for _, _, num in fig_refs], [84, 85, 86, 87])

    def test_figures_reports_permissions_and_exercise_numerical_alignment(self):
        self.assertEqual(verify_ch00_figures_and_reports(ROOT), [])
        ex_data = ex00.run_experiments()["exercises"]

        for num, png_name, json_name in CH00_FIG_SPECS:
            png_mode = stat.S_IMODE((ROOT / "figures" / png_name).stat().st_mode)
            json_mode = stat.S_IMODE((ROOT / "codes/chapters/ch00/reports" / json_name).stat().st_mode)
            self.assertEqual(png_mode, 0o644)
            self.assertEqual(json_mode, 0o644)

        r84 = strict_json_loads((ROOT / "codes/chapters/ch00/reports/figure84_complex_phasor_fourier.json").read_bytes())
        self.assertAlmostEqual(r84["anchors"]["parseval_energy"], ex_data["E00-04"]["time_energy"], places=15)
        r85 = strict_json_loads((ROOT / "codes/chapters/ch00/reports/figure85_subspace_projection_svd.json").read_bytes())
        self.assertAlmostEqual(r85["anchors"]["condition_number_2"], ex_data["E00-07"]["condition_number_2"], places=12)
        r86 = strict_json_loads((ROOT / "codes/chapters/ch00/reports/figure86_hermitian_rayleigh_gev.json").read_bytes())
        self.assertAlmostEqual(r86["anchors"]["max_gev_sinr"], ex_data["E00-09"]["max_gev_sinr"], places=14)
        r87 = strict_json_loads((ROOT / "codes/chapters/ch00/reports/figure87_wirtinger_probability_bayesian.json").read_bytes())
        self.assertAlmostEqual(r87["anchors"]["mvdr_output_power"], ex_data["E00-11"]["mvdr_output_power"], places=14)

    def test_negative_controls_detect_tampered_types_digests_and_sources(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_root = Path(tmpdir)
            (tmp_root / "figures").mkdir(parents=True)
            (tmp_root / "codes/chapters/ch00/reports").mkdir(parents=True)
            for rel in EXPECTED_CH00_PLOT_SOURCES:
                dst = tmp_root / rel
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(ROOT / rel, dst)
            for _, png_name, json_name in CH00_FIG_SPECS:
                shutil.copy2(ROOT / "figures" / png_name, tmp_root / "figures" / png_name)
                shutil.copy2(
                    ROOT / "codes/chapters/ch00/reports" / json_name,
                    tmp_root / "codes/chapters/ch00/reports" / json_name,
                )
            self.assertEqual(verify_ch00_figures_and_reports(tmp_root), [])

            # Tamper 1: float figure number 84.0 or bool anchor value rejected even if PNG digest were updated
            rep84_path = tmp_root / "codes/chapters/ch00/reports/figure84_complex_phasor_fourier.json"
            orig84 = json.loads(rep84_path.read_text(encoding="utf-8"))
            bad_fig = copy.deepcopy(orig84)
            bad_fig["figure"] = 84.0
            rep84_path.write_text(json.dumps(bad_fig), encoding="utf-8")
            self.assertTrue(verify_ch00_figures_and_reports(tmp_root))

            bad_bool = copy.deepcopy(orig84)
            bad_bool["anchors"]["parseval_energy"] = True
            rep84_path.write_text(json.dumps(bad_bool), encoding="utf-8")
            self.assertTrue(verify_ch00_figures_and_reports(tmp_root))

            # Restore report, tamper source file bytes -> caught by source_sha256 check
            shutil.copy2(
                ROOT / "codes/chapters/ch00/reports/figure84_complex_phasor_fourier.json",
                rep84_path,
            )
            for rel in ("codes/chapters/ch00/core/chapter00_math.py",
                        "codes/chapters/ch02/core/conventions.py"):
                with self.subTest(source=rel):
                    path = tmp_root / rel
                    path.write_text("# modified\n", encoding="utf-8")
                    errs = verify_ch00_figures_and_reports(tmp_root)
                    self.assertTrue(any("source SHA mismatch" in e for e in errs))
                    shutil.copy2(ROOT / rel, path)


if __name__ == "__main__":
    unittest.main()
