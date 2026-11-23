"""Captured HTML bytes, dated rebuilds and parsed PDF provenance controls."""
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from pypdf import PdfReader, PdfWriter

from scripts import build_pdf
from scripts.pdf_input_identity import (PDFInputIdentity, html_input_identity,
    require_unchanged_pdf_input, verify_canonical_html_input, verify_pdf_input_metadata)


SOURCE = '0123456789ab'


def html(date='2026-11-23', body='Original formula and prose'):
    return (f'<html><head><meta name="tutorial-source-digest" content="{SOURCE}">'
            f'<meta name="tutorial-build-date" content="{date}"></head>'
            f'<body>{body}</body></html>').encode()


def pdf(identity):
    writer = PdfWriter()
    writer.add_blank_page(200, 200)
    writer.add_metadata(identity.metadata())
    output = io.BytesIO()
    writer.write(output)
    return output.getvalue()


class PDFInputIdentityTests(unittest.TestCase):
    def test_final_validation_cannot_publish_a_changed_input(self):
        for change in ('private', 'public', 'source'):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as directory:
                out = Path(directory)
                combined = out / 'combined.html'
                combined.write_bytes(html())
                target = out / 'microphone-array-tutorial.pdf'
                previous = pdf(html_input_identity(html()))
                target.write_bytes(previous)
                state = {'digest': SOURCE}

                def render(private, destination, *args, input_identity):
                    state['private'] = private
                    destination.write_bytes(pdf(input_identity))

                def final_links(destination):
                    if change == 'source':
                        state['digest'] = '111111111111'
                    else:
                        changed = combined if change == 'public' else state['private']
                        changed.write_bytes(html(body='Changed during final validation'))

                patches = {
                    'OUT': out, 'source_digest': lambda: state['digest'],
                    'build_html': lambda date: (html(date).decode(), []),
                    'outline_from_html': lambda value: [],
                    'check_figures': lambda: None, 'check_mathjax_assets': lambda: None,
                    'print_pdf': render, 'validate_pdf_links': final_links,
                }
                from contextlib import ExitStack
                with ExitStack() as stack:
                    for name, value in patches.items():
                        stack.enter_context(patch.object(build_pdf, name, value))
                    args = ['--no-bookmarks', '--min-pages', '1']
                    args += ['--pdf-only'] if change == 'public' else ['--build-date', '2026-11-23']
                    with self.assertRaises(SystemExit):
                        build_pdf.main(args)
                self.assertEqual(target.read_bytes(), previous)

    def test_body_and_script_metadata_do_not_replace_the_single_head_identity(self):
        value = html(body='<meta name="tutorial-source-digest" content="bad">'
                          '<script>"<meta name=\"tutorial-build-date\" content=\"bad\">"</script>')
        self.assertEqual(html_input_identity(value).source_digest, SOURCE)
        for invalid in (value.replace(b'</head>', b'<meta name="tutorial-source-digest" content="bad"></head>'),
                        value.replace(b'2026-11-23', b'2026-13-23'),
                        value.replace(b'name="tutorial-source-digest"', b'name="tutorial-source-digest" name="other"'),
                        value.replace(b'<head>', b'').replace(b'</head>', b'')):
            with self.assertRaises(ValueError):
                html_input_identity(invalid)

    def test_valid_source_marker_cannot_hide_an_edited_body(self):
        with self.assertRaisesRegex(ValueError, 'canonical body'):
            verify_canonical_html_input(html(body='Edited formula'), lambda date: html(date).decode(), lambda: SOURCE)
        identity = verify_canonical_html_input(html(), lambda date: (html(date).decode(), []), lambda: SOURCE)
        self.assertEqual(identity, html_input_identity(html()))

    def test_pdf_build_date_controls_a_later_validation(self):
        identity = html_input_identity(html('2020-01-02'))
        dates = []
        def builder(date):
            dates.append(date)
            return html(date).decode(), []
        result = verify_pdf_input_metadata(PdfReader(io.BytesIO(pdf(identity))), builder, lambda: SOURCE)
        self.assertEqual(result, identity)
        self.assertEqual(dates, ['2020-01-02'])

    def test_wrong_input_sha_and_wrong_metadata_types_are_rejected(self):
        identity = html_input_identity(html())
        changed = PDFInputIdentity(SOURCE, '0' * 64, identity.build_date)
        with self.assertRaisesRegex(ValueError, 'canonical dated HTML'):
            verify_pdf_input_metadata(PdfReader(io.BytesIO(pdf(changed))), lambda date: html(date).decode(), lambda: SOURCE)
        for metadata in ({}, {'/SourceDigest': SOURCE, '/HTMLInputSHA256': 3, '/BuildDate': '2026-11-23'}):
            with self.assertRaises(ValueError):
                PDFInputIdentity.from_metadata(metadata)

    def test_finite_source_and_private_snapshot_change_guards(self):
        values = iter((SOURCE, '111111111111'))
        with self.assertRaisesRegex(ValueError, 'changed during'):
            verify_canonical_html_input(html(), lambda date: html(date).decode(), lambda: next(values))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'private.html'
            path.write_bytes(html())
            identity = html_input_identity(html())
            require_unchanged_pdf_input(identity, path, lambda: SOURCE)
            path.write_bytes(html(body='Changed private input'))
            with self.assertRaisesRegex(ValueError, 'private PDF HTML'):
                require_unchanged_pdf_input(identity, path, lambda: SOURCE)

    def test_bookmark_arguments_cannot_relabel_an_existing_pdf(self):
        original = html_input_identity(html(body='Old source body'))
        replacement = html_input_identity(html())
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'printed.pdf'
            path.write_bytes(pdf(original))
            before = path.read_bytes()
            with patch.object(build_pdf, 'source_digest', return_value=SOURCE):
                with self.assertRaisesRegex(ValueError, 'another HTML input'):
                    build_pdf.add_bookmarks(path, [], input_identity=replacement)
            self.assertEqual(path.read_bytes(), before)


if __name__ == '__main__':
    unittest.main()
