"""Regression controls for selected English publication text at full size."""
import unittest
from pathlib import Path
import tempfile
from unittest.mock import patch

from pypdf import PdfReader, PdfWriter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from scripts import make_figures_en as figures
from scripts import build_pdf, build_pdf_en as builder
from scripts import build_site_en


class ReviewedFigureFontTests(unittest.TestCase):
    def tearDown(self):
        plt.close('all')

    def test_selected_legend_keeps_eleven_points(self):
        fig,axes=plt.subplots(2,2,figsize=(9.5,8.5),constrained_layout=True)
        ax=axes[1,0]
        for text in ('Confidence','Ungated mask','Fused mask','Lip occlusion'):
            ax.plot([0,1],[0,1],label=text)
        ax.legend(ncol=2)
        protected=figures._publication_font_artists(fig,'fig83_av_asynchrony_occlusion.png')
        figures._autofit_translated_layout(fig,set(),protected)
        figures._validate_publication_fonts(fig,'fig83_av_asynchrony_occlusion.png',protected)
        self.assertTrue(all(t.get_fontsize()==11 for t in ax.get_legend().get_texts()))

    def test_long_reviewed_title_rejects_instead_of_shrinking(self):
        fig,axes=plt.subplots(3,1,figsize=(9.5,10))
        axes[2].set_title('A very long translated label '*30)
        protected=figures._publication_font_artists(fig,'fig12_gcc_phat.png')
        figures._autofit_translated_layout(fig,set(),protected)
        self.assertEqual(axes[2].title.get_fontsize(),12)
        with self.assertRaisesRegex(ValueError,'shorten or wrap'):
            figures._validate_publication_fonts(fig,'fig12_gcc_phat.png',protected)

    def test_unrelated_figure_keeps_existing_fit_policy(self):
        fig,ax=plt.subplots()
        ax.set_title('An unrelated label '*20,fontsize=12)
        protected=figures._publication_font_artists(fig,'fig01_binaural.png')
        self.assertEqual(protected,{})
        figures._autofit_translated_layout(fig,set(),protected)
        self.assertLess(ax.title.get_fontsize(),12)

    def test_regenerated_report_cannot_borrow_old_canonical_digest(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);canonical=root/'canonical.json';staged=root/'staged.json'
            canonical.write_bytes(b'{"value":1}\n');staged.write_bytes(canonical.read_bytes())
            self.assertEqual(figures._staged_report_bytes(canonical,staged,Path.read_bytes),staged.read_bytes())
            staged.write_bytes(b'{"value":2}\n')
            with self.assertRaisesRegex(ValueError,'rebuild Chinese figures first'):
                figures._staged_report_bytes(canonical,staged,Path.read_bytes)
            canonical.unlink()
            with self.assertRaises(ValueError):figures._staged_report_bytes(canonical,staged,Path.read_bytes)


class EnglishPDFInputIdentityTests(unittest.TestCase):
    DIGEST='a'*12
    DATE='2001-02-03'

    def html(self,date):
        return ('<!doctype html><html><head>'
            f'<meta name="tutorial-source-digest" content="{self.DIGEST}">'
            f'<meta name="tutorial-build-date" content="{date}">'
            '</head><body>Canonical formula body</body></html>')

    def run_pdf_only(self,folder,print_callback):
        with patch.object(builder,'OUT',folder),patch.object(builder,'source_digest_en',return_value=self.DIGEST),\
             patch.object(builder,'build_html_en',side_effect=lambda date:(self.html(date),[])),\
             patch.object(builder,'print_pdf_en',side_effect=print_callback),\
             patch.object(build_pdf,'check_mathjax_assets'),patch.object(build_pdf,'outline_from_html',return_value=[]),\
             patch.object(build_pdf,'validate_pdf_links'),\
             patch('scripts.quality_check_en.cached_html_issues_en',return_value=[]):
            builder.main(['--pdf-only','--no-bookmarks','--min-pages','0'])

    def test_old_date_pdf_uses_private_snapshot_and_no_bookmarks_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            folder=Path(directory);cached=folder/'combined_en.html';cached.write_text(self.html(self.DATE))
            captured={}
            def print_pdf(snapshot,destination,_minimum,*,input_identity):
                self.assertNotEqual(snapshot,cached)
                self.assertEqual(snapshot.read_bytes(),self.html(self.DATE).encode())
                self.assertEqual(input_identity.build_date,self.DATE)
                captured['identity']=input_identity
                writer=PdfWriter();writer.add_blank_page(width=100,height=100)
                writer.add_metadata(input_identity.metadata());writer.write(destination)
                return 1
            self.run_pdf_only(folder,print_pdf)
            reader=PdfReader(folder/'microphone-array-tutorial-en.pdf')
            self.assertEqual(build_pdf.PDFInputIdentity.from_metadata(reader.metadata),captured['identity'])
            self.assertEqual(build_pdf.verify_pdf_input_metadata(reader,lambda date:self.html(date),lambda:self.DIGEST),captured['identity'])

    def test_cached_body_edit_is_rejected_even_with_valid_meta(self):
        with tempfile.TemporaryDirectory() as directory:
            folder=Path(directory);(folder/'combined_en.html').write_text(self.html(self.DATE).replace('formula','altered'))
            called=[]
            with self.assertRaisesRegex(SystemExit,'canonical body'):
                self.run_pdf_only(folder,lambda *a,**k:called.append(True))
            self.assertEqual(called,[])

    def test_private_snapshot_mutation_prevents_publication(self):
        with tempfile.TemporaryDirectory() as directory:
            folder=Path(directory);(folder/'combined_en.html').write_text(self.html(self.DATE))
            def print_pdf(snapshot,destination,_minimum,*,input_identity):
                writer=PdfWriter();writer.add_blank_page(width=100,height=100)
                writer.add_metadata(input_identity.metadata());writer.write(destination)
                snapshot.write_text('changed private input')
                return 1
            with self.assertRaisesRegex(SystemExit,'private PDF HTML input changed'):
                self.run_pdf_only(folder,print_pdf)
            self.assertFalse((folder/'microphone-array-tutorial-en.pdf').exists())

    def test_public_cache_change_is_preserved_and_prevents_publication(self):
        with tempfile.TemporaryDirectory() as directory:
            folder=Path(directory);cached=folder/'combined_en.html';cached.write_text(self.html(self.DATE))
            def print_pdf(snapshot,destination,_minimum,*,input_identity):
                self.assertEqual(snapshot.read_text(),self.html(self.DATE))
                writer=PdfWriter();writer.add_blank_page(width=100,height=100)
                writer.add_metadata(input_identity.metadata());writer.write(destination)
                cached.write_text(self.html(self.DATE).replace('Canonical','Changed'))
                return 1
            with self.assertRaises(SystemExit):self.run_pdf_only(folder,print_pdf)
            self.assertIn('Changed',cached.read_text())
            self.assertFalse((folder/'microphone-array-tutorial-en.pdf').exists())

    def test_final_guard_follows_heavy_metadata_validation(self):
        original_verify=build_pdf.verify_pdf_input_metadata
        for kind in ('private HTML','source'):
            with self.subTest(kind=kind),tempfile.TemporaryDirectory() as directory:
                folder=Path(directory);(folder/'combined_en.html').write_text(self.html(self.DATE))
                capture={}
                def print_pdf(snapshot,destination,_minimum,*,input_identity):
                    capture['snapshot']=snapshot
                    writer=PdfWriter();writer.add_blank_page(width=100,height=100)
                    writer.add_metadata(input_identity.metadata());writer.write(destination)
                    return 1
                def verify(reader,html_builder,digest):
                    result=original_verify(reader,html_builder,digest)
                    if kind=='private HTML':capture['snapshot'].write_text('changed after heavy validation')
                    else:digest.return_value='b'*12
                    return result
                with patch.object(build_pdf,'verify_pdf_input_metadata',side_effect=verify),self.assertRaises(SystemExit):
                    self.run_pdf_only(folder,print_pdf)
                self.assertFalse((folder/'microphone-array-tutorial-en.pdf').exists())

    def test_bookmark_caller_cannot_resign_an_old_pdf(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'input.pdf'
            old=build_pdf.html_input_identity(self.html(self.DATE).encode())
            new=build_pdf.PDFInputIdentity(old.source_digest,'b'*64,old.build_date)
            writer=PdfWriter();writer.add_blank_page(width=100,height=100)
            writer.add_metadata(old.metadata());writer.write(path);before=path.read_bytes()
            with patch.object(builder,'source_digest_en',return_value=self.DIGEST),\
                 patch.object(builder,'PdfWriter',side_effect=AssertionError('must reject before cloning')),\
                 self.assertRaisesRegex(ValueError,'bookmark input identity differs'):
                builder.add_bookmarks_en(path,[],input_identity=new)
            self.assertEqual(path.read_bytes(),before)


class CompleteEnglishDependencyTests(unittest.TestCase):
    def test_all_report_bytes_and_real_shared_oracle_sources_are_read(self):
        from scripts.figure_report_contracts import FIGURE_REPORTS
        original=Path.read_bytes
        common=['scripts/figure_report_contracts.py','scripts/portable_figure_report_checks.py',
            'scripts/quality_check.py','scripts/pdf_input_identity.py']
        common += [entry.report_relative for entry in FIGURE_REPORTS]
        root=Path(__file__).resolve().parents[1]
        for digest in (build_site_en.source_digest_en,builder.source_digest_en):
            before=digest()
            for relative in common:
                reads=[];target=root/relative
                def changed(path):
                    data=original(path)
                    if path==target:reads.append(path);data+=b'independent dependency edit'
                    return data
                with self.subTest(builder=digest.__module__,dependency=relative),patch.object(Path,'read_bytes',changed):
                    self.assertNotEqual(digest(),before)
                self.assertEqual(reads,[target])
            self.assertEqual(digest(),before)

    def test_printed_room_result_image_bytes_bind_english_pdf(self):
        root=Path(__file__).resolve().parents[1]
        target=root/'codes/chapters/appendix_b/room_audio/ROOM_RESULTS.png'
        before=builder.source_digest_en();original=Path.read_bytes;reads=[]
        def changed(path):
            data=original(path)
            if path==target:reads.append(path);data+=b'independent image edit'
            return data
        with patch.object(Path,'read_bytes',changed):self.assertNotEqual(builder.source_digest_en(),before)
        self.assertEqual(reads,[target]);self.assertEqual(builder.source_digest_en(),before)


if __name__=='__main__':unittest.main()
