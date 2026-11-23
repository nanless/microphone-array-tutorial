"""Independent rejection controls for English publication source identity."""
import hashlib
import unittest
from pathlib import Path
from unittest.mock import patch

from scripts import build_pdf_en, build_site_en, publication_dependencies_en as deps

ROOT = Path(__file__).resolve().parents[1]


class EnglishDependencyTests(unittest.TestCase):
    def test_consumed_dependency_changes_invalidate_both_publishers(self):
        original = Path.read_bytes
        common = (
            "scripts/chapter_identity.py", "scripts/chapter_numbering.json",
            "scripts/legacy_sequential_anchors.json", "scripts/heading_aliases.py",
            "scripts/build_site.py", "scripts/build_markdown_helpers.py",
            "scripts/inline_layout.js", "scripts/make_figures_en.py",
            "scripts/figure_en_translations.json", "chapters/00_math-foundations.md",
            "scripts/make_aec_figures.py", "scripts/make_wearable_figures.py",
            "codes/chapters/ch00/research/01_spatial_and_tracking.md",
        )
        cases = ((build_site_en.source_digest_en, common + (
            "scripts/course_home_en/home.js", "assets/readme/cases/locata/MANIFEST.json",
            "assets/readme/cases/locata/NOTICE.txt", "assets/readme/cases/locata/before.wav",
            "assets/readme/lab/array-lab.html", "codes/chapters/ch05/phase_audio/phase_flip.wav")),
            (build_pdf_en.source_digest_en, common + (
                "scripts/build_pdf.py", "scripts/pdf_math.js",
                "scripts/vendor/mathjax-3.2.2/tex-svg.js")))
        for digest, paths in cases:
            before = digest()
            self.assertEqual(before, digest())
            for relative in paths:
                target = ROOT / relative
                def changed(path):
                    return original(path) + (b"independent edit" if path == target else b"")
                with self.subTest(publisher=digest.__module__, dependency=relative):
                    with patch.object(Path, "read_bytes", changed):
                        self.assertNotEqual(before, digest())

    def test_pdf_only_vendor_and_unconsumed_file_boundaries(self):
        site = deps.site_sources(ROOT, build_site_en.SRC_EN, build_site_en.RESEARCH_EN,
                                 [name for name, _ in build_site_en.RESEARCH_EN_PAGES])
        pdf = deps.pdf_sources(ROOT, build_pdf_en.SRC_EN, build_site_en.RESEARCH_EN,
                               [name for name, _ in build_site_en.RESEARCH_EN_PAGES])
        vendor = ROOT / "scripts/vendor/mathjax-3.2.2/tex-svg.js"
        self.assertNotIn(vendor, site)
        self.assertIn(vendor, pdf)
        self.assertNotIn(ROOT / "reviews/2026-11-23-multiround-review.md", site + pdf)

    def test_digest_binds_path_and_bytes_and_ignores_input_order(self):
        # Independent in-memory oracle, not the production path enumerator.
        a, b = ROOT / "a", ROOT / "b"
        contents = {a: b"AB", b: b"C"}
        expected = hashlib.sha256(b"a\0AB\0b\0C\0").hexdigest()[:12]
        with patch.object(Path, "read_bytes", lambda p: contents[p]):
            self.assertEqual(deps.source_digest(ROOT, [b, a, a]), expected)
            self.assertNotEqual(deps.source_digest(ROOT, [a]), expected)



class EnglishOrphanSpaceTests(unittest.TestCase):
    BLOCK = b'/Span << /ActualText ( ) >> BDC 1 0 0 1 20 80 Tm (\\000\\003) Tj EMC'

    def test_inherited_font_must_be_proved_and_visible_glyph_stays_content(self):
        from tests.test_pdf_layout_spaces import page_fixture, operations, plain_operations
        for visible, count in ((False, 1), (True, 0)):
            writer, page, _tree, _font = page_fixture(
                b'BT /Fblank 20 Tf '+self.BLOCK+b' ET', visible_space=visible)
            before = operations(page, writer)
            self.assertEqual(build_pdf_en._mark_orphan_layout_span_spaces(writer), count)
            after = operations(page, writer)
            self.assertEqual(plain_operations([x for x in before if x[1] != b'BDC']),
                             plain_operations([x for x in after if x[1] != b'BDC']))
            if not count:
                self.assertEqual(plain_operations(before), plain_operations(after))

    def test_unknown_font_wrong_unicode_and_zero_advance_stay_semantic(self):
        from tests.test_pdf_layout_spaces import page_fixture, unicode_map
        cases = ({}, {'cmap': unicode_map(b'<0003> <0041>', count=1)},
                 {'space_advance': 0}, {'visible_space': True})
        for options in cases:
            with self.subTest(options=options):
                writer, page, _tree, _font = page_fixture(
                    b'BT /Fblank 20 Tf '+self.BLOCK+b' ET', **options)
                if not options:
                    del page['/Resources']['/Font']['/Fblank']
                self.assertEqual(build_pdf_en._mark_orphan_layout_span_spaces(writer), 0)

    def test_q_Q_restore_actual_font_and_protect_semantic_marks(self):
        from pypdf.generic import NameObject
        from tests.test_pdf_layout_spaces import page_fixture, font_resource, operations
        writer, page, _tree, _font = page_fixture(
            b'BT /Fblank 20 Tf q /Fvisible 20 Tf '+self.BLOCK+b' Q '+self.BLOCK+b' ET')
        page['/Resources']['/Font'][NameObject('/Fvisible')] = font_resource(writer, visible_space=True)
        self.assertEqual(build_pdf_en._mark_orphan_layout_span_spaces(writer), 1)
        tags = [op[0][0] for op in operations(page, writer) if op[1] == b'BDC']
        self.assertEqual(tags, ['/Span', '/Artifact'])
        for outer in (b'/P << /MCID 0 >> BDC ', b'/Artifact BMC '):
            writer, page, _tree, _font = page_fixture(
                b'BT /Fblank 20 Tf '+outer+self.BLOCK+b' EMC ET')
            self.assertEqual(build_pdf_en._mark_orphan_layout_span_spaces(writer), 0)

    def test_strings_that_look_like_pdf_markers_are_not_parsed_as_operators(self):
        from tests.test_pdf_layout_spaces import page_fixture
        writer, _page, _tree, _font = page_fixture(
            b'BT /Fblank 20 Tf ( /Span << /ActualText \\( \\) >> BDC EMC ) Tj ET')
        self.assertEqual(build_pdf_en._mark_orphan_layout_span_spaces(writer), 0)

    def test_semantic_change_preserves_unrelated_high_precision_coordinates(self):
        from tests.test_pdf_layout_spaces import page_fixture
        raw = (b'BT /Fblank 20 Tf ' + self.BLOCK + b' ET\nq\n'
               b'1000000000 0 0 1000000000 -1234567800 -1234567800 cm\n'
               b'1.2345678912345 1.2345678912345 0.00000005 0.00000005 re f\nQ\n')
        writer, page, _tree, _font = page_fixture(raw)
        self.assertEqual(build_pdf_en._mark_orphan_layout_span_spaces(writer), 1)
        actual = page.get_contents().get_data()
        self.assertIn(b'1.2345678912345 1.2345678912345', actual)
        self.assertIn(b'0.00000005 0.00000005 re', actual)

class EnglishPngTests(unittest.TestCase):
    def test_decoded_identity_translation_and_wrong_same_generator_figure(self):
        import tempfile
        from PIL import Image, PngImagePlugin
        from scripts import quality_check_en as q
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            zh, en = root/'zh', root/'en'
            zh.mkdir(); en.mkdir()
            name = 'fig01_binaural.png'
            metadata = {
                'SourceScript': 'scripts/make_figures.py',
                'SourceScriptDigest': hashlib.sha256((ROOT/'scripts/make_figures.py').read_bytes()).hexdigest(),
                'FigureFilename': name,
                'TranslationScript': 'scripts/make_figures_en.py',
                'TranslationScriptSHA256': hashlib.sha256((ROOT/'scripts/make_figures_en.py').read_bytes()).hexdigest(),
                'TranslationTableSHA256': hashlib.sha256((ROOT/'scripts/figure_en_translations.json').read_bytes()).hexdigest(),
            }
            def save(path, data):
                info = PngImagePlugin.PngInfo()
                for k,v in data.items(): info.add_text(k,v)
                Image.new('RGB', (800,300), 'white').save(path, pnginfo=info)
            save(zh/name, metadata)
            save(en/name, metadata)
            with patch.object(q, 'FIGURES_ZH', zh):
                self.assertEqual(q.english_png_issues(en/name, name), [])
                for field in ('FigureFilename','SourceScriptDigest','TranslationScript',
                              'TranslationScriptSHA256','TranslationTableSHA256'):
                    bad = dict(metadata); bad[field] = 'fig02_dsin_geometry.png' if field=='FigureFilename' else 'stale'
                    save(en/name, bad)
                    with self.subTest(field=field):
                        self.assertTrue(any(field in x for x in q.english_png_issues(en/name, name)))
                (en/name).write_bytes(b'\x89PNG\r\n\x1a\n'+bytes(1100))
                self.assertTrue(q.english_png_issues(en/name, name))
                save(en/name, metadata)
                raw = (en/name).read_bytes()
                # Valid PNG header and metadata, but incomplete pixel stream.
                offset = raw.index(b'IDAT')
                (en/name).write_bytes(raw[:offset+5])
                self.assertTrue(q.english_png_issues(en/name, name))

    def test_actual_english_saver_adds_identity_and_retains_numeric_metadata(self):
        from scripts import make_figures_en as m
        expected = m.translation_metadata('fig55_mask_representation.png')
        self.assertEqual(expected['FigureFilename'], 'fig55_mask_representation.png')
        self.assertEqual(expected['TranslationScript'], 'scripts/make_figures_en.py')
        self.assertEqual(expected['TranslationTableSHA256'], hashlib.sha256(m.TRANSLATIONS_PATH.read_bytes()).hexdigest())
        import tempfile
        from PIL import Image
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory)
            def draw(*_args, **_kwargs):
                figure, axis = m.plt.subplots()
                axis.plot([0, 1], [0, 1])
                figure.savefig(m.FIGURES_ZH_DIR/'fig55_mask_representation.png', metadata={
                    'SourceScript': 'scripts/make_figures.py',
                    'NumericalReportDigest': 'independent report identity',
                    'AudioManifestDigest': 'independent audio identity'})
                m.plt.close(figure)
            with patch.object(m, 'FIGURES_EN_DIR', target), patch.object(m.runpy, 'run_path', side_effect=draw):
                written = m.generate_all_english_figures(('make_figures.py',))
            self.assertEqual(written, [target/'fig55_mask_representation.png'])
            with Image.open(written[0]) as picture:
                picture.load()
                self.assertEqual(picture.info['FigureFilename'], 'fig55_mask_representation.png')
                self.assertEqual(picture.info['TranslationTableSHA256'], expected['TranslationTableSHA256'])
                self.assertEqual(picture.info['TranslationScriptSHA256'], hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest())
                self.assertEqual(picture.info['NumericalReportDigest'], 'independent report identity')
                self.assertEqual(picture.info['AudioManifestDigest'], 'independent audio identity')



class EnglishOutlineTests(unittest.TestCase):
    EXPECTED = [('Chapter One', 'ch-1', [('Repeated title', 'ch-1-sec-1-1', [])])]

    def test_cached_html_checks_digest_independent_ids_and_duplicates(self):
        from scripts import quality_check_en as q
        html = ('<html><head><meta name="tutorial-build-date" content="2001-02-03">'
                '<meta name="tutorial-source-digest" content="0123456789ab"></head>'
                '<body>Source sha256 0123456789ab'
                '<div class="chap" id="ch-1"><h1>Chapter One</h1>'
                '<h2 id="ch-1-sec-1-1">Repeated title</h2><p>body</p></div></body>')
        with patch.object(q, 'expected_pdf_outline_en', return_value=self.EXPECTED):
            self.assertEqual(q.cached_html_issues_en(html, '0123456789ab'), [])
            for bad in (html.replace('0123456789ab','fedcba987654'),
                        html.replace('Source sha256 0123456789ab',''),
                        html.replace('ch-1-sec-1-1','ch-1-sec-1-2'),
                        html.replace('<p>body</p>','<p id="ch-1-sec-1-1">body</p>')):
                self.assertTrue(q.cached_html_issues_en(bad, '0123456789ab'))

    def test_pdf_outlines_follow_named_XYZ_even_with_repeated_titles(self):
        import io
        from pypdf import PdfReader, PdfWriter
        from pypdf.generic import ArrayObject, FloatObject, NameObject, TextStringObject, Fit
        from scripts import quality_check_en as q
        def reader(*, wrong_page=False, wrong_top=False, empty=False):
            writer = PdfWriter(); pages=[writer.add_blank_page(200,200) for _ in range(2)]
            for name,top in [('ch-1',180),('ch-1-sec-1-1',150)]:
                writer.add_named_destination_array(TextStringObject(name),ArrayObject([
                    pages[1].indirect_reference,NameObject('/XYZ'),FloatObject(20),FloatObject(top),FloatObject(1)]))
            if not empty:
                parent=writer.add_outline_item('Chapter One',1,fit=Fit.xyz(20,180,1))
                writer.add_outline_item('Repeated title',0 if wrong_page else 1,parent=parent,
                                        fit=Fit.xyz(20,140 if wrong_top else 150,1))
            data=io.BytesIO();writer.write(data);return PdfReader(io.BytesIO(data.getvalue()))
        with patch.object(q,'expected_pdf_outline_en',return_value=self.EXPECTED):
            self.assertEqual(q.pdf_outline_issues_en(reader()),[])
            for options in ({'wrong_page':True},{'wrong_top':True},{'empty':True}):
                with self.subTest(options=options):self.assertTrue(q.pdf_outline_issues_en(reader(**options)))

    def test_expectation_comes_from_source_and_independent_catalog(self):
        from scripts import quality_check_en as q
        catalog=(('01_fixture.md','Independent Chapter'),)
        documents={'01_fixture.md':'## Ignored source title\n### 1.1 A section\n#### A detail\n'}
        with patch.object(q,'EXPECTED_PDF_CHAPTERS_EN',catalog), \
             patch.object(q.quality_check,'EXPECTED_SECTION_COUNTS',{'01_fixture.md':1}), \
             patch.object(q.quality_check,'EXPECTED_SUBSECTION_COUNTS',{'01_fixture.md':1}):
            value=q.expected_pdf_outline_en(documents)
            self.assertEqual(value[0][0:2],('Independent Chapter','ch-1'))
            self.assertEqual(value[0][2][0][0:2],('1.1 A section','ch-1-sec-1-1'))
            with self.assertRaises(ValueError):q.expected_pdf_outline_en({'01_fixture.md':'## No section\n'})


class EnglishPdfPreflightTests(unittest.TestCase):
    def test_pdf_only_rejects_stale_missing_or_wrong_content_before_printer(self):
        import tempfile
        from scripts import quality_check_en as q
        html = ('<html><head><meta name="tutorial-build-date" content="2001-02-03">'
                '<meta name="tutorial-source-digest" content="0123456789ab"></head>'
                '<body>Source sha256 0123456789ab'
                '<div class="chap" id="ch-1"><h1>Chapter One</h1>'
                '<h2 id="ch-1-sec-1-1">Repeated title</h2><p>body</p></div></body>')
        cases = ((html, True), (html.replace('0123456789ab','fedcba987654'),False),
                 (html.replace('Source sha256 0123456789ab',''),False),
                 (html.replace('ch-1-sec-1-1','ch-1-wrong'),False))
        for content, reaches_printer in cases:
            with self.subTest(reaches_printer=reaches_printer), tempfile.TemporaryDirectory() as directory:
                out=Path(directory);(out/'combined_en.html').write_text(content)
                with patch.object(build_pdf_en,'OUT',out), \
                     patch.object(build_pdf_en,'source_digest_en',return_value='0123456789ab'), \
                     patch.object(build_pdf_en,'build_html_en',return_value=(html,EnglishOutlineTests.EXPECTED)), \
                     patch.object(q,'expected_pdf_outline_en',return_value=EnglishOutlineTests.EXPECTED), \
                     patch.object(build_pdf_en.build_pdf,'check_mathjax_assets'), \
                     patch.object(build_pdf_en,'print_pdf_en',side_effect=RuntimeError('printer sentinel')) as printer:
                    with self.assertRaises(RuntimeError if reaches_printer else (SystemExit,ValueError)):
                        build_pdf_en.main(['--pdf-only'])
                    self.assertEqual(printer.call_count, int(reaches_printer))

    def test_pdf_gate_rejects_actual_reader_metadata_even_with_other_controls_unchanged(self):
        import io,tempfile
        from pypdf import PdfReader,PdfWriter
        from scripts import quality_check_en as q
        writer=PdfWriter();writer.add_blank_page(200,200)
        writer.add_metadata({'/SourceDigest':'0123456789ab'})
        output=io.BytesIO();writer.write(output)
        with tempfile.TemporaryDirectory() as directory:
            dist=Path(directory)
            for name in ('microphone-array-tutorial.pdf','microphone-array-tutorial-en.pdf'):
                (dist/name).write_bytes(output.getvalue())
            def check(digest):
                def reader(_path):
                    value=PdfReader(io.BytesIO(output.getvalue()))
                    from pypdf.generic import NameObject, TextStringObject
                    if digest is None:del value.trailer['/Info']['/SourceDigest']
                    else:value.trailer['/Info'][NameObject('/SourceDigest')]=TextStringObject(digest)
                    return value
                errors=[]
                with patch.object(q,'DIST',dist),patch.object(q,'PdfReader',side_effect=reader), \
                     patch.object(q,'pdf_outline_issues_en',return_value=[]), \
                     patch.object(build_pdf_en,'source_digest_en',return_value='0123456789ab'), \
                     patch.object(q.build_pdf,'source_digest',return_value='0123456789ab'), \
                     patch.object(q.build_pdf,'validate_pdf_structure_tags'), \
                     patch.object(q.build_pdf,'validate_pdf_links'):
                    q.check_pdfs(errors)
                # Small untagged test PDFs intentionally fail unrelated page/end gates.
                return [error for error in errors if 'SourceDigest' in error]
            self.assertEqual(check('0123456789ab'),[])
            self.assertEqual(len(check('fedcba987654')),2)
            self.assertEqual(len(check(None)),2)


if __name__ == "__main__":
    unittest.main()
