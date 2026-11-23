"""PNG export uses the chapter report root independently of its output layout."""
import hashlib
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from PIL import Image

from scripts import make_figures as drawing
from scripts.figure_report_contracts import REPORT_BY_NUMBER
from scripts.make_figures_en import _staged_report_bytes


class FigureSaveContractTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name) / 'workspace'
        self.root.mkdir()

    def figure(self):
        fig, ax = drawing.plt.subplots(figsize=(1.5, 1.5))
        ax.plot([0, 1], [0, 1])
        self.addCleanup(drawing.plt.close, fig)
        return fig

    def report(self, chapters, number=75, payload=b'{"value": 1}\n'):
        path = chapters / Path(REPORT_BY_NUMBER[number].report_relative).relative_to('codes/chapters')
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
        return path

    def save(self, out, chapters, name, metadata=None):
        out.mkdir(parents=True, exist_ok=True)
        with patch.object(drawing, 'OUT', out), patch.object(drawing, 'CODE_CHAPTERS', chapters):
            drawing.save(self.figure(), name, extra_metadata=metadata)
        return out / name

    def metadata(self, path):
        with Image.open(path) as image:
            image.load()
            return image.info.copy()

    def test_generic_and_unregistered_numbered_png_keep_source_metadata(self):
        for name in ('metadata-smoke.png', 'fig01_ordinary.png'):
            with self.subTest(name=name):
                path = self.save(self.root / 'out', self.root / 'absent-chapters', name)
                info = self.metadata(path)
                self.assertEqual(info['SourceScript'], 'scripts/make_figures.py')
                self.assertEqual(info['SourceScriptDigest'], drawing.source_script_digest())
                self.assertEqual(info['FigureFilename'], name)
                self.assertNotIn('NumericalReportDigest', info)

    def test_report_root_is_independent_of_output_parent_and_layout(self):
        for out in (self.root, self.root / 'figures', self.root / 'nested' / 'output'):
            with self.subTest(out=out):
                chapters = self.root / 'authoritative' / 'codes' / 'chapters'
                report = self.report(chapters)
                # The old parent-derived location contains a different report.
                self.report(out.parent / 'codes' / 'chapters', payload=b'{"decoy": true}\n')
                info = self.metadata(self.save(out, chapters, REPORT_BY_NUMBER[75].png_filename))
                self.assertEqual(info['NumericalReportDigest'], hashlib.sha256(report.read_bytes()).hexdigest())

    def test_default_relative_layout_and_current_report_bytes(self):
        chapters = self.root / 'codes' / 'chapters'
        out = self.root / 'figures'
        report = self.report(chapters)
        name = REPORT_BY_NUMBER[75].png_filename
        first = self.metadata(self.save(out, chapters, name))['NumericalReportDigest']
        report.write_bytes(b'{"value": 2}\n')
        second = self.metadata(self.save(out, chapters, name))['NumericalReportDigest']
        self.assertNotEqual(first, second)
        self.assertEqual(second, hashlib.sha256(report.read_bytes()).hexdigest())

    def test_registered_filename_variants_are_rejected_before_png_write(self):
        chapters = self.root / 'codes' / 'chapters'
        self.report(chapters)
        for name in ('fig75_wrong.png', 'fig75.png', 'fig75wrong.png', 'fig075_imaging_objectives.png'):
            with self.subTest(name=name):
                with self.assertRaisesRegex(ValueError, 'filename'):
                    self.save(self.root / 'out', chapters, name)
                self.assertFalse((self.root / 'out' / name).exists())

    def test_provided_digest_must_match_real_report_and_cannot_replace_it(self):
        chapters = self.root / 'codes' / 'chapters'
        report = self.report(chapters)
        name = REPORT_BY_NUMBER[75].png_filename
        expected = hashlib.sha256(report.read_bytes()).hexdigest()
        path = self.save(self.root / 'out', chapters, name, {'NumericalReportDigest': expected})
        self.assertEqual(self.metadata(path)['NumericalReportDigest'], expected)
        path.unlink()
        with self.assertRaisesRegex(ValueError, 'digest'):
            self.save(self.root / 'out', chapters, name, {'NumericalReportDigest': '0' * 64})
        self.assertFalse(path.exists())
        report.unlink()
        with self.assertRaises(FileNotFoundError):
            self.save(self.root / 'out', chapters, name, {'NumericalReportDigest': expected})
        self.assertFalse(path.exists())

    def test_english_staged_read_binds_actual_equal_bytes_and_rejects_difference(self):
        chapters = self.root / 'canonical' / 'codes' / 'chapters'
        canonical = self.report(chapters)
        staged = self.report(self.root / 'staged' / 'codes' / 'chapters')
        original_read = Path.read_bytes

        def staged_read(path):
            if path == canonical:
                return _staged_report_bytes(canonical, staged, original_read)
            return original_read(path)

        name = REPORT_BY_NUMBER[75].png_filename
        with patch.object(Path, 'read_bytes', staged_read):
            path = self.save(self.root / 'english', chapters, name)
        self.assertEqual(self.metadata(path)['NumericalReportDigest'], hashlib.sha256(original_read(staged)).hexdigest())
        path.unlink()
        staged.write_bytes(b'{"different": true}\n')
        with patch.object(Path, 'read_bytes', staged_read):
            with self.assertRaisesRegex(ValueError, 'regenerated report differs'):
                self.save(self.root / 'english', chapters, name)
        self.assertFalse(path.exists())


if __name__ == '__main__':
    unittest.main()
