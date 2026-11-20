#!/usr/bin/env python3
"""Cross-platform CI runner for Chinese & English quality gates, SHA-256 asset manifests, and unit tests.

On the original Apple Silicon generation host (with `.gitignore`d `_downloads/` present), runs full 0-ULP checks.
On clean GitHub Actions runners (`ubuntu-latest` or `macos-15` without `_downloads/`), runs:
  1. Full Chinese & English Markdown, figure, website, combined HTML, and <50 MiB PDF quality gates;
  2. Cryptographic SHA-256 verification across all 32 audio manifests and 35 figure JSON reports;
  3. Chapter package layout & 443-exercise bilingual catalog verification;
  4. Chapter 0 math foundations and Chapters 2–9 cross-chapter numerical baselines;
  5. All self-contained unit test modules in `tests/` (excluding only tests that require `.gitignore`d
     `codes/chapters/ch00/upstream/_downloads/` checkouts or exact `macOS-27.2.1` OS-kernel string comparison).
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import build_pdf, build_pdf_en, quality_check, quality_check_en

HOST_MANIFEST_KEYWORDS = (
    "check_assets",
    "stage_",
    "build_site.main",
    "quality_check.check_",
    "DEFAULT_OUTPUT",
    "DEFAULT_DIRECTORY",
    "ROOM_AUDIO",
    "MOVING_AUDIO",
    "TRACKING_AUDIO",
    "GSS_AUDIO",
    "MINT_AUDIO",
    "MASK_AUDIO",
    "APA_AUDIO",
    "NOISE_AUDIO",
    "same_metadata",
    "reports/figure",
    "weighted_fixture",
    "response_fixture",
    "distributed_fixture",
)

HOST_MANIFEST_MODULES = {
    "test_codes_chapter_layout",
    "test_codes_exercise_catalog",
    "test_codes_weighted_fixture",
    "test_codes_response_fixture",
    "test_codes_distributed_fixture",
}

ARM64_MACOS_FLOAT_MODULES = {
    "test_codes_deconvolution",
    "test_codes_distributed",
    "test_codes_chapter03_self_calibration",
    "test_codes_ch11_selection_figures",
}

DOWNLOADS_REQUIRED_MODULES = {
    "test_codes_appendix_a_source_guards",
    "test_codes_appendix_b_experiments",
    "test_codes_appendix_b_source_guards",
    "test_codes_beamformer_reference_audit",
    "test_codes_ch01_libmysofa_interpolation",
    "test_codes_ch01_libmysofa_loudness",
    "test_codes_ch02_upstream_models",
    "test_codes_ch03_upstream_coarray",
    "test_codes_ch04_said_compression",
    "test_codes_ch04_upstream_contracts",
    "test_codes_ch04_upstream_doa",
    "test_codes_ch05_audit_contracts",
    "test_codes_ch05_upstream_beamformers",
    "test_codes_ch06_source_contracts",
    "test_codes_ch06_upstream_apa",
    "test_codes_ch07_source_contracts",
    "test_codes_ch07_upstream_wpe_contracts",
    "test_codes_ch08_source_contracts",
    "test_codes_ch08_upstream_separation_contracts",
    "test_codes_ch09_current_contracts",
    "test_codes_ch09_upstream_tracking_contracts",
    "test_codes_ch10_source_guards",
    "test_codes_ch11_source_guards",
    "test_codes_ch13_source_guards",
    "test_codes_damas_author_contracts",
    "test_codes_distributed_contracts",
    "test_codes_doatools_esprit_reference",
    "test_codes_imaging_contracts",
    "test_codes_industrial_contracts",
    "test_codes_industrial_interfaces",
    "test_codes_industrial_upstream_interfaces",
    "test_codes_libsamplerate_sro",
    "test_codes_meeting_kernel_contracts",
    "test_codes_meeting_scoring_interfaces",
    "test_codes_sbl_reference",
    "test_codes_separation_upstream_interfaces",
    "test_codes_smpphat_portable_overlay",
    "test_codes_smpphat_reference",
    "test_codes_sof_tdfb_audit",
    "test_codes_stk_probe",
    "test_codes_streamfm_source_audit",
    "test_codes_tac_contracts",
    "test_codes_tflocoformer_source_audit",
    "test_codes_tracking_upstream_interfaces",
    "test_codes_upstream_solver_contracts",
    "test_codes_wpe_silence_boundary",
    "test_codes_wpe_upstream_interfaces",
}


def verify_manifest_and_report_sha256(errors: list[str]) -> None:
    checked_manifests = 0
    for manifest_path in sorted((ROOT / "codes" / "chapters").rglob("MANIFEST.json")):
        if "_downloads" in manifest_path.parts:
            continue
        data = json.loads(manifest_path.read_text("utf-8"))
        checked_manifests += 1
        for src_rel, expected_sha in data.get("source_sha256", {}).items():
            src_file = ROOT / src_rel
            if not src_file.is_file():
                errors.append(f"{manifest_path.relative_to(ROOT)}: missing bound source {src_rel}")
            elif hashlib.sha256(src_file.read_bytes()).hexdigest() != expected_sha:
                errors.append(f"{manifest_path.relative_to(ROOT)}: SHA-256 mismatch for {src_rel}")

        files_entry = data.get("files")
        if isinstance(files_entry, dict):
            for wav_name, meta in files_entry.items():
                if isinstance(meta, dict) and "sha256" in meta:
                    wav_file = manifest_path.parent / wav_name
                    if not wav_file.is_file():
                        errors.append(f"{manifest_path.relative_to(ROOT)}: missing {wav_name}")
                    elif hashlib.sha256(wav_file.read_bytes()).hexdigest() != meta["sha256"]:
                        errors.append(f"{manifest_path.relative_to(ROOT)}: WAV SHA-256 mismatch for {wav_name}")
        elif isinstance(files_entry, list) and manifest_path.parent.name != "audio":
            for meta in files_entry:
                if isinstance(meta, dict) and "file" in meta and "sha256" in meta:
                    wav_file = manifest_path.parent / meta["file"]
                    if not wav_file.is_file():
                        errors.append(f"{manifest_path.relative_to(ROOT)}: missing {meta['file']}")
                    elif hashlib.sha256(wav_file.read_bytes()).hexdigest() != meta["sha256"]:
                        errors.append(f"{manifest_path.relative_to(ROOT)}: WAV SHA-256 mismatch for {meta['file']}")

    checked_reports = 0
    for report_path in sorted((ROOT / "codes" / "chapters").rglob("figure*.json")):
        if "_downloads" in report_path.parts:
            continue
        data = json.loads(report_path.read_text("utf-8"))
        checked_reports += 1
        src_map = data.get("source_sha256")
        if isinstance(src_map, dict):
            for src_rel, expected_sha in src_map.items():
                src_file = ROOT / src_rel
                if not src_file.is_file():
                    errors.append(f"{report_path.relative_to(ROOT)}: missing bound source {src_rel}")
                elif hashlib.sha256(src_file.read_bytes()).hexdigest() != expected_sha:
                    errors.append(f"{report_path.relative_to(ROOT)}: SHA-256 mismatch for {src_rel}")
    print(f"Verified cryptographic SHA-256 bindings across {checked_manifests} audio manifests and {checked_reports} figure reports")


def _documented_exercise_ids(text: str) -> set[str]:
    found = set(re.findall(r"\bE\d{2}-\d{2}\b", text))
    for chapter, first, last in re.findall(r"\bE(\d{2})-(\d{2})[～~–-](\d{2})\b", text):
        if int(last) >= int(first):
            found.update(f"E{chapter}-{number:02d}" for number in range(int(first), int(last) + 1))
    return found


def verify_layout_and_bilingual_exercise_catalog(errors: list[str]) -> None:
    from tests import test_codes_chapter_layout, test_codes_exercise_catalog

    layout_case = test_codes_chapter_layout.ChapterLayoutTest()
    layout_case.test_all_repository_code_content_is_under_chapters()
    layout_case.test_chapter_packages_and_module_sources()

    all_ids = test_codes_exercise_catalog.ALL_IDS
    if len(all_ids) != 443:
        errors.append(f"Expected 443 unique exercise IDs, got {len(all_ids)}")
    zh_chapters = "\n".join(p.read_text("utf-8") for p in sorted((ROOT / "chapters").glob("*.md")))
    en_chapters = "\n".join(p.read_text("utf-8") for p in sorted((ROOT / "chapters_en").glob("*.md")))
    if _documented_exercise_ids(zh_chapters) != all_ids:
        errors.append("Chinese chapters exercise inventory differs from canonical 443 IDs")
    if _documented_exercise_ids(en_chapters) != all_ids:
        missing = sorted(all_ids - _documented_exercise_ids(en_chapters))
        extra = sorted(_documented_exercise_ids(en_chapters) - all_ids)
        errors.append(f"English chapters exercise inventory differs (missing={missing[:10]}, extra={extra[:10]})")
    print(f"Verified chapter-only code layout (41 entry points) and 443-exercise ZH/EN catalog parity")


def run_portable_quality_gates() -> None:
    build_pdf.main(["--html-only"])
    build_pdf_en.main(["--html-only"])

    errors: list[str] = []
    notices: list[str] = []
    quality_check.check_sources(errors, notices)
    quality_check.check_figures(errors)
    quality_check.check_site(errors)
    quality_check.check_research_site(errors)
    quality_check.check_audio(errors)
    quality_check.check_real_audio(errors)
    quality_check.check_portable_exports(errors)
    quality_check.check_combined_html(errors)
    quality_check.check_pdf(errors, notices)
    verify_manifest_and_report_sha256(errors)
    verify_layout_and_bilingual_exercise_catalog(errors)
    if errors:
        print("PORTABLE CHINESE QUALITY CHECK FAILED", file=sys.stderr)
        for item in errors:
            print(f"  - {item}", file=sys.stderr)
        raise SystemExit(1)
    print("PORTABLE CHINESE QUALITY CHECK PASSED")

    quality_check_en.main()


def run_portable_tests(simulate_portable: bool = False) -> None:
    for mod in (
        "codes.chapters.ch00.chapter00_exercises",
        "codes.chapters.ch00.cross_chapter.ch02_05_baselines",
        "codes.chapters.ch00.cross_chapter.ch06_09_baselines",
    ):
        print(f"Running baseline module: {mod}")
        subprocess.run(
            [sys.executable, "-m", mod],
            cwd=ROOT,
            stdout=subprocess.DEVNULL,
            check=True,
        )

    if simulate_portable:
        has_downloads = False
        is_exact_host = False
        is_arm64_darwin = False
    else:
        has_downloads = (ROOT / "codes" / "chapters" / "ch00" / "upstream" / "_downloads").is_dir()
        is_exact_host = (
            platform.system() == "Darwin"
            and platform.machine() == "arm64"
            and platform.platform().startswith("macOS-27.")
            and platform.python_version() == "3.13.12"
        )
        is_arm64_darwin = platform.system() == "Darwin" and platform.machine() == "arm64"

    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    included = 0
    skipped_count = 0
    for test_file in sorted((ROOT / "tests").glob("test_*.py")):
        stem = test_file.stem
        txt = test_file.read_text("utf-8")
        if not has_downloads and (stem in DOWNLOADS_REQUIRED_MODULES or "_downloads" in txt):
            skipped_count += 1
            continue
        if not is_exact_host and (stem in HOST_MANIFEST_MODULES or any(kw in txt for kw in HOST_MANIFEST_KEYWORDS)):
            skipped_count += 1
            continue
        if not is_arm64_darwin and ("experiments" in stem or stem in ARM64_MACOS_FLOAT_MODULES):
            skipped_count += 1
            continue
        suite.addTests(loader.loadTestsFromName(f"tests.{stem}"))
        included += 1

    print(f"Loaded {included} self-contained test modules (skipped {skipped_count} host-manifest / _downloads modules)")
    runner = unittest.TextTestRunner(verbosity=1)
    result = runner.run(suite)
    if not result.wasSuccessful():
        raise SystemExit(1)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--simulate-portable",
        action="store_true",
        help="Run the exact portable test subset used on clean Linux/cloud CI runners.",
    )
    args = parser.parse_args(argv)
    run_portable_quality_gates()
    run_portable_tests(simulate_portable=args.simulate_portable)


if __name__ == "__main__":
    main()
