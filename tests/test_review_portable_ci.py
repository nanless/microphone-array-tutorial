"""Independent capability decisions and the verified Pages staging path."""
import builtins
import os
from pathlib import Path, PureWindowsPath
import subprocess
import sys
import unittest
from unittest.mock import patch

from scripts import run_portable_ci as ci


class PortableSelectionTests(unittest.TestCase):
    PORTABLE = dict(has_downloads=False, is_exact_host=False, is_arm64_darwin=False)

    def test_each_reason_is_visible_and_multiple_reasons_are_retained(self):
        self.assertEqual(ci.test_selection_reasons('test_plain_math', 'assert 1+1 == 2', **self.PORTABLE), [])
        self.assertTrue(any('experiment filename' in r for r in ci.test_selection_reasons(
            'test_plain_experiments', '', **self.PORTABLE)))
        reasons = ci.test_selection_reasons('test_codes_appendix_b_experiments',
            '# _downloads and check_assets', **self.PORTABLE)
        self.assertEqual(len(reasons), 4)
        self.assertTrue(any('conservative text' in r for r in reasons))
        self.assertTrue(any('generation-host source' in r for r in reasons))

    def test_original_host_and_intermediate_capabilities(self):
        self.assertEqual(ci.test_selection_reasons('test_codes_appendix_b_experiments',
            '_downloads check_assets', has_downloads=True, is_exact_host=True, is_arm64_darwin=True), [])
        reasons = ci.test_selection_reasons('test_plain_experiments', 'check_assets',
            has_downloads=True, is_exact_host=False, is_arm64_darwin=True)
        self.assertEqual(len(reasons), 1)
        self.assertIn('generation-host source', reasons[0])

    def test_policy_fixtures_are_selected_without_enabling_similar_real_tests(self):
        text = '_downloads check_assets'
        self.assertEqual(ci.test_selection_reasons('test_review_portable_ci', text, **self.PORTABLE), [])
        self.assertTrue(ci.test_selection_reasons('test_real_portable_ci', text, **self.PORTABLE))

    def test_pages_checks_before_staging_in_the_same_job(self):
        source = (Path(__file__).resolve().parents[1]/'.github/workflows/pages.yml').read_text()
        build, deploy = source.split('\n  deploy:', 1)
        self.assertIn('run: python scripts/run_portable_ci.py', build)
        self.assertLess(build.index('run: python scripts/run_portable_ci.py'), build.index('mkdir -p _pages'))
        self.assertLess(build.index('mkdir -p _pages'), build.index('actions/upload-pages-artifact'))
        self.assertEqual(build.count('actions/checkout@'), 1)
        self.assertIn('needs: build', deploy)
        self.assertIn("python-version: '3.13'", build)
        self.assertIn("node-version: '20'", build)

    def test_redirected_console_can_report_chinese_in_legacy_locale(self):
        root = Path(__file__).resolve().parents[1]
        probe = """
import runpy
from unittest.mock import patch
from scripts import build_pdf

def announce(_args):
    print("合订：19 篇")
    raise SystemExit(0)

with patch.object(build_pdf, "main", side_effect=announce):
    runpy.run_path("scripts/run_portable_ci.py", run_name="__main__")
"""
        result = subprocess.run(
            [sys.executable, "-B", "-c", probe], cwd=root,
            env={**os.environ, "PYTHONIOENCODING": "cp1252:strict"},
            capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr.decode("utf-8", errors="replace"))
        self.assertEqual(result.stdout.decode("utf-8").splitlines(), ["合订：19 篇"])


    def test_publication_identity_is_independent_of_windows_path_order(self):
        from scripts import build_pdf, build_pdf_en, build_site, build_site_en
        from scripts import publication_dependencies_en as dependencies

        calls = (build_site.source_digest, build_pdf.source_digest,
                 build_site_en.source_digest_en, build_pdf_en.source_digest_en)
        platform_order = []

        def windows_sorted(values, *, key=None, reverse=False):
            paths = list(values)
            if key is None and paths and all(isinstance(path, Path) for path in paths):
                key = lambda path: PureWindowsPath(path.as_posix())
                platform_order.append(len(paths))
            return builtins.sorted(paths, key=key, reverse=reverse)

        expected = tuple(check() for check in calls)
        with (
            patch.object(build_site, 'sorted', windows_sorted, create=True),
            patch.object(build_pdf, 'sorted', windows_sorted, create=True),
            patch.object(dependencies, 'sorted', windows_sorted, create=True),
        ):
            current = tuple(check() for check in calls)
        self.assertTrue(platform_order)
        self.assertEqual(current, expected)


if __name__ == '__main__':
    unittest.main()
