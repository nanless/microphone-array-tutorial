#!/usr/bin/env python3
"""Portable publication identities, 35 numeric reports, and selected unit tests.

Check bilingual sources/sites/PDFs, all32 audio manifests, and fixed35 report
identities, independent arithmetic and report-to-PNG bindings. Small analytic
and actual-PCM integer checks use explicit tolerances where appropriate; saved
random events are checked as observations and summaries, not resimulated trials.
This layer does not promise full waveform 0-ULP replay on Linux or other hosts.

The test subset depends on checkout availability, original generation-host
identity and Apple Silicon arithmetic requirements. Additional strict replay
tests are selected on the original host. Every excluded module and reason is
printed; their absence must not be described as successful original execution.
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
import wave

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import build_pdf, build_pdf_en, quality_check, quality_check_en
from codes.chapters.ch00.io_contracts import strict_json_loads, validate_parent_chain

# These reviewed fixtures use capability flags, private HTML/PDF examples or
# independent report/PCM arithmetic, without upstream imports or full waveform
# replay. Their source intentionally names capabilities and asset constants.
PORTABLE_POLICY_FIXTURES = {
    "test_review_portable_ci", "test_review_r4_portable_reports",
    "test_review_r4_publication_en",
}

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

# These exercise tests isolate strict published-asset replay per test method.
EXPERIMENTS_WITH_LOCAL_AUDIO_FIXTURES = {
    "test_codes_chapter01_experiments",
    "test_codes_chapter02_experiments",
    "test_codes_chapter03_experiments",
    "test_codes_chapter04_experiments",
    "test_codes_chapter05_experiments",
    "test_codes_chapter07_experiments",
    "test_codes_chapter08_experiments",
    "test_codes_chapter11_experiments",
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


EXPECTED_AUDIO_MANIFESTS = (
    'codes/chapters/appendix_a/weighted_audio/MANIFEST.json',
    'codes/chapters/appendix_b/response_audio/MANIFEST.json',
    'codes/chapters/appendix_b/room_audio/MANIFEST.json',
    'codes/chapters/ch00/audio/MANIFEST.json',
    'codes/chapters/ch01/binaural_audio/MANIFEST.json',
    'codes/chapters/ch01/spectral_audio/MANIFEST.json',
    'codes/chapters/ch02/real_audio/MANIFEST.json',
    'codes/chapters/ch02/stft_audio/MANIFEST.json',
    'codes/chapters/ch02/sweep_audio/MANIFEST.json',
    'codes/chapters/ch03/baseline_audio/MANIFEST.json',
    'codes/chapters/ch03/geometry_audio/MANIFEST.json',
    'codes/chapters/ch04/focus_audio/MANIFEST.json',
    'codes/chapters/ch04/reflection_audio/MANIFEST.json',
    'codes/chapters/ch05/binaural_beamforming_audio/MANIFEST.json',
    'codes/chapters/ch05/derivative_audio/MANIFEST.json',
    'codes/chapters/ch05/phase_audio/MANIFEST.json',
    'codes/chapters/ch06/apa_audio/MANIFEST.json',
    'codes/chapters/ch06/reference_audio/MANIFEST.json',
    'codes/chapters/ch07/delay_audio/MANIFEST.json',
    'codes/chapters/ch07/mint_audio/MANIFEST.json',
    'codes/chapters/ch08/css_audio/MANIFEST.json',
    'codes/chapters/ch08/gss_audio/MANIFEST.json',
    'codes/chapters/ch08/mask_audio/MANIFEST.json',
    'codes/chapters/ch09/moving_audio/MANIFEST.json',
    'codes/chapters/ch09/tracking_audio/MANIFEST.json',
    'codes/chapters/ch10/channel_audio/MANIFEST.json',
    'codes/chapters/ch10/noise_audio/MANIFEST.json',
    'codes/chapters/ch11/scenario_audio/MANIFEST.json',
    'codes/chapters/ch12/imaging_audio/MANIFEST.json',
    'codes/chapters/ch13/distributed_audio/MANIFEST.json',
    'codes/chapters/ch14/wearable_audio/MANIFEST.json',
    'codes/chapters/ch15/tse_audio/MANIFEST.json',
)

def _check_bound_file(folder, filename, digest, label, errors):
    try:
        if (type(filename) is not str or not filename or filename in (".", "..")
                or "/" in filename or "\\" in filename):
            raise ValueError("bound member must be a simple filename")
        if type(digest) is not str or re.fullmatch(r"[0-9a-f]{64}", digest) is None:
            raise ValueError("bound SHA-256 must be64 lowercase hexadecimal characters")
        target = validate_parent_chain(folder / filename)
        if not target.is_file() or hashlib.sha256(target.read_bytes()).hexdigest() != digest:
            raise ValueError("actual member SHA-256 differs")
    except (OSError, ValueError, TypeError) as exc:
        errors.append(f"{label}/{filename}: {exc}")


def check_binaural_measurement_bindings(folder, data, errors):
    """Check this manifest's four measurement members, without full synthesis replay."""
    expected = {role: "binaural_" + role + ".wav" for role in ("input", "bmvdr", "blcmv", "null")}
    measurements = data.get("measurements")
    if type(measurements) is not dict or set(measurements) != set(expected):
        errors.append("binaural measurement member roles differ")
        return
    for role, filename in expected.items():
        row = measurements[role]
        if type(row) is not dict or row.get("file") != filename:
            errors.append("binaural measurement filename differs: " + role)
            continue
        _check_bound_file(folder, filename, row.get("sha256"), "binaural measurements", errors)
        try:
            with wave.open(str(validate_parent_chain(folder / filename)), "rb") as wav:
                if (wav.getframerate(), wav.getnchannels(), wav.getsampwidth(),
                        wav.getnframes(), wav.getcomptype()) != (16000, 2, 2, 32000, "NONE"):
                    raise ValueError("binaural PCM format differs")
                if len(wav.readframes(32000)) != 128000:
                    raise ValueError("binaural PCM payload is truncated")
        except (OSError, ValueError, EOFError, wave.Error) as exc:
            errors.append(f"{filename}: {exc}")


def verify_manifest_and_report_sha256(errors: list[str]) -> None:
    from scripts.figure_report_contracts import FIGURE_REPORTS
    declared = {ROOT / name for name in EXPECTED_AUDIO_MANIFESTS}
    actual = {p for p in (ROOT / "codes/chapters").rglob("MANIFEST.json") if "_downloads" not in p.parts}
    if actual != declared:
        errors.append("audio manifest inventory differs from the fixed32 publication members")
    checked_manifests = 0
    for manifest_path in sorted(declared):
        label = manifest_path.relative_to(ROOT).as_posix()
        try:
            data = strict_json_loads(validate_parent_chain(manifest_path).read_bytes())
            if type(data) is not dict:
                raise ValueError("manifest root must be an object")
            checked_manifests += 1
            sources = data.get("source_sha256", {})
            if type(sources) is not dict:
                raise ValueError("source_sha256 must be an object")
            for relative, expected_sha in sources.items():
                relative_path = Path(relative)
                if relative_path.is_absolute() or ".." in relative_path.parts or "\\" in relative:
                    raise ValueError("bound source must stay within the repository")
                source = validate_parent_chain(ROOT / relative)
                if type(expected_sha) is not str or re.fullmatch(r"[0-9a-f]{64}", expected_sha) is None:
                    raise ValueError("invalid bound source digest")
                if not source.is_file() or hashlib.sha256(source.read_bytes()).hexdigest() != expected_sha:
                    errors.append(f"{label}: source SHA-256 mismatch for {relative}")
            files = data.get("files")
            if isinstance(files, dict):
                for filename, meta in files.items():
                    if not isinstance(meta, dict) or "sha256" not in meta:
                        raise ValueError("bound file record is missing its SHA-256")
                    _check_bound_file(manifest_path.parent, filename, meta["sha256"], label, errors)
            elif isinstance(files, list) and manifest_path.parent.name != "audio":
                for meta in files:
                    if not isinstance(meta, dict) or "file" not in meta or "sha256" not in meta:
                        raise ValueError("bound file list record is missing filename/SHA-256")
                    _check_bound_file(manifest_path.parent, meta["file"], meta["sha256"], label, errors)
            elif label != "codes/chapters/ch05/binaural_beamforming_audio/MANIFEST.json" and manifest_path.parent.name != "audio":
                raise ValueError("manifest files must be a dictionary or list")
            if label == "codes/chapters/ch05/binaural_beamforming_audio/MANIFEST.json":
                check_binaural_measurement_bindings(manifest_path.parent, data, errors)
        except (OSError, ValueError, TypeError, KeyError) as exc:
            errors.append(f"{label}: {exc}")
    checked_reports = 0
    for contract in FIGURE_REPORTS:
        path = ROOT / contract.report_relative
        try:
            data = strict_json_loads(validate_parent_chain(path).read_bytes())
            if type(data) is not dict:
                raise ValueError("report root must be an object")
            checked_reports += 1
            sources = data.get("source_sha256")
            if isinstance(sources, dict):
                for relative, expected_sha in sources.items():
                    candidate = Path(relative)
                    if candidate.is_absolute() or ".." in candidate.parts or "\\" in relative:
                        raise ValueError("report source must stay within the repository")
                    if hashlib.sha256(validate_parent_chain(ROOT / relative).read_bytes()).hexdigest() != expected_sha:
                        errors.append(f"{contract.report_relative}: source SHA-256 mismatch for {relative}")
        except (OSError, ValueError, TypeError, KeyError) as exc:
            errors.append(f"{contract.report_relative}: {exc}")
    print(f"Parsed {checked_manifests}/32 audio manifests and {checked_reports}/35 figure reports; "
          "checked declared dictionary SHA bindings and actual member bytes. "
          "Independent35 report schema/source-role/arithmetic checks run separately.")




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
    if len(all_ids) != 360:
        errors.append(f"Expected 360 unique exercise IDs, got {len(all_ids)}")
    zh_chapters = "\n".join(p.read_text("utf-8") for p in sorted((ROOT / "chapters").glob("*.md")))
    en_chapters = "\n".join(p.read_text("utf-8") for p in sorted((ROOT / "chapters_en").glob("*.md")))
    if _documented_exercise_ids(zh_chapters) != all_ids:
        errors.append("Chinese chapters exercise inventory differs from canonical 360 IDs")
    if _documented_exercise_ids(en_chapters) != all_ids:
        missing = sorted(all_ids - _documented_exercise_ids(en_chapters))
        extra = sorted(_documented_exercise_ids(en_chapters) - all_ids)
        errors.append(f"English chapters exercise inventory differs (missing={missing[:10]}, extra={extra[:10]})")
    print(f"Verified chapter-only code layout (41 entry points) and 360-exercise ZH/EN catalog parity")


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
    from scripts.portable_figure_report_checks import check_figure_reports
    check_figure_reports(errors)
    verify_layout_and_bilingual_exercise_catalog(errors)
    if errors:
        print("PORTABLE CHINESE QUALITY CHECK FAILED", file=sys.stderr)
        for item in errors:
            print(f"  - {item}", file=sys.stderr)
        raise SystemExit(1)
    print("PORTABLE CHINESE QUALITY CHECK PASSED")

    quality_check_en.main()


def is_original_generation_host() -> bool:
    return (platform.system() == "Darwin"
            and platform.machine() == "arm64"
            and platform.platform().startswith("macOS-27.")
            and platform.python_version() == "3.13.12")


def test_selection_reasons(stem: str, text: str, *, has_downloads: bool,
                           is_exact_host: bool, is_arm64_darwin: bool) -> list[str]:
    """Return all matching conservative selection rules without importing a test."""
    if stem in PORTABLE_POLICY_FIXTURES:
        return []
    reasons = []
    if not has_downloads:
        if stem in DOWNLOADS_REQUIRED_MODULES:
            reasons.append("requires local upstream checkouts (explicit module)")
        if "_downloads" in text:
            reasons.append("source mentions _downloads (conservative text rule)")
    if not is_exact_host:
        if stem in HOST_MANIFEST_MODULES:
            reasons.append("requires original generation host (explicit module)")
        if stem not in EXPERIMENTS_WITH_LOCAL_AUDIO_FIXTURES:
            matches = [word for word in HOST_MANIFEST_KEYWORDS if word in text]
            if matches:
                reasons.append("generation-host source rule: " + ", ".join(matches))
    if not is_arm64_darwin:
        if "experiments" in stem:
            reasons.append("experiment filename requires Apple Silicon macOS")
        if stem in ARM64_MACOS_FLOAT_MODULES:
            reasons.append("Apple Silicon macOS floating-point baseline (explicit module)")
    return reasons


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
        is_exact_host = is_original_generation_host()
        is_arm64_darwin = platform.system() == "Darwin" and platform.machine() == "arm64"

    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    included = 0
    skipped_count = 0
    for test_file in sorted((ROOT / "tests").glob("test_*.py")):
        stem = test_file.stem
        txt = test_file.read_text("utf-8")
        reasons = test_selection_reasons(
            stem, txt, has_downloads=has_downloads, is_exact_host=is_exact_host,
            is_arm64_darwin=is_arm64_darwin)
        if reasons:
            skipped_count += 1
            print(f"Excluded tests.{stem}: {'; '.join(reasons)}")
            continue
        suite.addTests(loader.loadTestsFromName(f"tests.{stem}"))
        included += 1

    print(f"Loaded {included} test modules; excluded {skipped_count} modules by the reasons above")
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
    # A redirected Windows console can otherwise reject Chinese progress messages.
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
    main()
