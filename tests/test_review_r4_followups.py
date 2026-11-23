"""Finite PHD, complex-storage, captured-PDF and full-size legend controls."""
from copy import deepcopy
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pypdf import PdfReader, PdfWriter

from scripts import build_pdf, build_pdf_en, make_figures_en
from scripts import portable_figure_report_checks as checks
from scripts.figure_report_contracts import REPORT_BY_NUMBER


ROOT = Path(__file__).resolve().parents[1]


def report(number):
    return json.loads((ROOT / REPORT_BY_NUMBER[number].report_relative).read_text())


def complex_records(value):
    if isinstance(value, dict):
        if set(value) == {'real', 'imag'}:
            yield value
        else:
            for item in value.values():
                yield from complex_records(item)
    elif isinstance(value, list):
        for item in value:
            yield from complex_records(item)


class PortableReportFollowups(unittest.TestCase):
    def test_phd_fixed_shape_preserves_legal_backgrounds(self):
        original = report(21)
        checks.validate_numeric_report(21, original)
        grid = original['phd_illustration']['grid_deg']
        gaussian = [math.exp(-((x+25)/4)**2/2)/(4*math.sqrt(2*math.pi))
                    + math.exp(-((x-40)/5)**2/2)/(5*math.sqrt(2*math.pi)) for x in grid]
        for background in ([0.]*400, [.00025]*400, [.00049*i/399 for i in range(400)]):
            curve = [x+y for x,y in zip(gaussian, background)]
            area = sum((grid[i+1]-grid[i])*(curve[i+1]+curve[i])/2 for i in range(399))
            changed = deepcopy(original)
            changed['phd_illustration']['intensity_targets_per_degree'] = [2*x/area for x in curve]
            checks.validate_numeric_report(21, changed)
        for curve in (original['phd_illustration']['intensity_targets_per_degree'][::-1], [2/180]*400):
            changed = deepcopy(original)
            changed['phd_illustration']['intensity_targets_per_degree'] = curve
            with self.assertRaisesRegex(ValueError, 'Gaussian/background'):
                checks.validate_numeric_report(21, changed)

    def test_all_complex_report_families_reject_coercion_and_broadcasting(self):
        # Other reports store real/imaginary pairs in ordinary fixed-shape arrays.
        for number in (62, 63, 66, 75, 77):
            original = report(number)
            checks._finite_tree(original)
            for index, value in enumerate(complex_records(original)):
                for bad in ('short', 'Boolean', 'string', 'ragged'):
                    changed = deepcopy(value)
                    if bad == 'short':
                        changed['imag'] = changed['imag'][:1]
                    elif bad == 'ragged':
                        changed['imag'] = [[0.], [0., 0.]]
                    else:
                        target = changed['imag']
                        while isinstance(target[0], list):
                            target = target[0]
                        target[0] = True if bad == 'Boolean' else '0.0'
                    with self.subTest(figure=number, record=index, bad=bad), self.assertRaises(ValueError):
                        checks._finite_tree(changed)

    def test_distributed_role_dimensions_and_startup_coordinates(self):
        original = report(66)
        checks._distributed_complex_shapes(original)
        for schedule in original['results']['trajectories']:
            solution = original['results']['trajectories'][schedule]['history'][0]['solutions'][0]
            for field in ('proposal_local_weights', 'new_compression', 'weights',
                          'effective_weights_after_broadcast', 'compressed_weights',
                          'raw_compressed_weights', 'normalized_projection', 'raw_projection'):
                changed = deepcopy(original)
                stored = changed['results']['trajectories'][schedule]['history'][0]['solutions'][0][field]
                for part in ('real', 'imag'):
                    stored[part] = solution[field][part][:1]
                with self.subTest(schedule=schedule, field=field), self.assertRaises(ValueError):
                    checks._distributed_complex_shapes(changed)
        changed = deepcopy(original)
        row = changed['results']['trajectories']['round_robin']['history'][0]
        for part in ('real', 'imag'):
            row['receiver_coefficients_raw_coordinates'][1][part].append(0.)
        with self.assertRaises(ValueError):
            checks._distributed_complex_shapes(changed)
        changed = deepcopy(original)
        del changed['results']['trajectories']['simultaneous']['history'][0]['solutions'][0]['incoming_compressions'][1]
        with self.assertRaises(ValueError):
            checks._distributed_complex_shapes(changed)


class CapturedEnglishPdfFollowups(unittest.TestCase):
    DIGEST = '0123456789ab'

    def html(self, date):
        return ('<html><head>'
                f'<meta name="tutorial-source-digest" content="{self.DIGEST}">'
                f'<meta name="tutorial-build-date" content="{date}"></head>'
                f'<body>Canonical dated body {date}</body></html>')

    def write_pdf(self, path, identity):
        writer = PdfWriter()
        writer.add_blank_page(100, 100)
        writer.add_metadata(identity.metadata())
        writer.write(path)

    def test_both_dates_match_capture_and_valid_old_pdf_cannot_replace_it(self):
        for date, swap in (('2001-02-03', False), ('2026-11-23', False), ('2026-11-23', True)):
            with self.subTest(date=date, swap=swap), tempfile.TemporaryDirectory() as directory:
                folder = Path(directory)
                output = folder/'microphone-array-tutorial-en.pdf'
                output.write_bytes(b'previous publication')
                captured = {}
                def printer(_html, target, _minimum, *, input_identity):
                    captured['identity'] = input_identity
                    self.write_pdf(target, input_identity)
                def links(target):
                    if swap:
                        self.write_pdf(target, build_pdf.html_input_identity(self.html('2001-02-03').encode()))
                with patch.object(build_pdf_en, 'OUT', folder), \
                     patch.object(build_pdf_en, 'source_digest_en', return_value=self.DIGEST), \
                     patch.object(build_pdf_en, 'build_html_en', side_effect=lambda day:(self.html(day), [])), \
                     patch.object(build_pdf_en, 'check_figures_en'), \
                     patch.object(build_pdf_en, 'print_pdf_en', side_effect=printer), \
                     patch.object(build_pdf, 'check_mathjax_assets'), \
                     patch.object(build_pdf, 'validate_pdf_links', side_effect=links):
                    if swap:
                        with self.assertRaisesRegex(SystemExit, 'captured HTML input'):
                            build_pdf_en.main(['--no-bookmarks', '--build-date', date, '--min-pages', '1'])
                        self.assertEqual(output.read_bytes(), b'previous publication')
                    else:
                        build_pdf_en.main(['--no-bookmarks', '--build-date', date, '--min-pages', '1'])
                        self.assertEqual(build_pdf.PDFInputIdentity.from_metadata(PdfReader(output).metadata), captured['identity'])


class Figure77LegendFollowups(unittest.TestCase):
    def tearDown(self):
        plt.close('all')

    def test_panel_b_legend_fits_at_eleven_points_and_overlong_text_rejects(self):
        translator = make_figures_en.FigureTranslator(make_figures_en.TRANSLATIONS_PATH)
        labels = [translator.translate_one(label) for label in
                  ('对侧头影区 210°～330°', '目标左耳幅度', 'BSM-LS', 'BSM-MagLS')]
        figure, axes = plt.subplots(2, 2, figsize=(9.5, 7.5), constrained_layout=True)
        axis = axes.flat[1]
        for label in labels:
            axis.plot([0, 1], [0, 1], label=label)
        axis.legend(ncol=2)
        protected = make_figures_en._publication_font_artists(figure, 'fig77_bsm_ls_vs_magls.png')
        make_figures_en._autofit_translated_layout(figure, set(), protected)
        make_figures_en._validate_publication_fonts(figure, 'fig77_bsm_ls_vs_magls.png', protected)
        self.assertEqual([text.get_fontsize() for text in axis.get_legend().get_texts()], [11.]*4)
        axis.get_legend().get_texts()[0].set_text('Overlong legend '*100)
        with self.assertRaisesRegex(ValueError, 'full-size legend'):
            make_figures_en._validate_publication_fonts(figure, 'fig77_bsm_ls_vs_magls.png', protected)
