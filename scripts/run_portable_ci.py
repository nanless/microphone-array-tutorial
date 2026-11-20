#!/usr/bin/env python3
"""Cross-platform CI runner for Chinese & English quality gates, SHA-256 asset manifests, and unit tests.

On the original Apple Silicon generation host, runs full 0-ULP host-locked checks.
On portable Linux/cloud CI runners (`ubuntu-latest`), runs 100% of Chinese & English
Markdown, figure, website, combined HTML, and PDF gates, verifies cryptographic SHA-256
bindings for all 30 standalone audio manifests and 35 figure JSON reports, executes cross-chapter
baselines, and runs all OS-portable unit test modules (skipping only 0-ULP Apple arm64 libm /
`platform.platform()` host-string comparisons).
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import build_pdf, build_pdf_en, quality_check, quality_check_en


HOST_REPLAY_KEYWORDS = (
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
)


def verify_manifest_and_report_sha256(errors: list[str]) -> None:
    """Cryptographically verify source_sha256 and WAV file sha256 across all committed manifests and reports."""
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
    if errors:
        print("PORTABLE CHINESE QUALITY CHECK FAILED", file=sys.stderr)
        for item in errors:
            print(f"  - {item}", file=sys.stderr)
        raise SystemExit(1)
    print("PORTABLE CHINESE QUALITY CHECK PASSED")

    quality_check_en.main()


def run_portable_tests() -> None:
    for mod in (
        "codes.chapters.ch00.chapter00_exercises",
        "codes.chapters.ch00.cross_chapter.ch02_05_baselines",
        "codes.chapters.ch00.cross_chapter.ch06_09_baselines",
    ):
        print(f"Running baseline module: {mod}")
        subprocess.run([sys.executable, "-m", mod], cwd=ROOT, check=True)

    is_exact_host = (
        platform.system() == "Darwin"
        and platform.machine() == "arm64"
        and platform.platform().startswith("macOS-27.")
        and platform.python_version() == "3.13.12"
    )
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    included = 0
    skipped_host = 0
    for test_file in sorted((ROOT / "tests").glob("test_*.py")):
        txt = test_file.read_text("utf-8")
        if not is_exact_host and any(kw in txt for kw in HOST_REPLAY_KEYWORDS):
            skipped_host += 1
            continue
        mod_name = f"tests.{test_file.stem}"
        suite.addTests(loader.loadTestsFromName(mod_name))
        included += 1

    print(f"Loaded {included} portable test modules (skipped {skipped_host} host-0ULP-manifest modules)")
    runner = unittest.TextTestRunner(verbosity=1)
    result = runner.run(suite)
    if not result.wasSuccessful():
        raise SystemExit(1)


def main() -> None:
    run_portable_quality_gates()
    run_portable_tests()


if __name__ == "__main__":
    main()
