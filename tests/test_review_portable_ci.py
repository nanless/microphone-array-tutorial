"""Independent capability decisions and the verified Pages staging path."""
import unittest
from pathlib import Path

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


if __name__ == '__main__':
    unittest.main()
