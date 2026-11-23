"""Authored PDF paths with extents independent of the vector optimizer."""

import io
import unittest

from pypdf import PdfReader, PdfWriter
from pypdf.generic import (ArrayObject, BooleanObject, ContentStream, DecodedStreamObject,
                          DictionaryObject, FloatObject, NameObject, NumberObject)

from scripts.build_pdf import (compress_pdf_page_streams_losslessly,
                               deduplicate_pdf_vector_glyphs, pdf_content_bytes_losslessly,
                               pdf_marked_true, preserve_pdf_float_values, tag_pdf_formulas,
                               verified_layout_space_font, write_pdf_with_object_streams)
from tests.test_pdf_layout_spaces import artifact_texts, font_resource, page_fixture


RECTANGLE = (b'80.000000 90.000000 m\n80.000000 90.000000 l\n'
             b'80.000000 90.000000 180.000000 150.000000 re\nh')
POLYGON = (b'100.000000 100.000000 m\n300.000000 100.000000 l\n'
           b'300.000000 300.000000 l\n100.000000 300.000000 l\nh')


def document(content):
    writer = PdfWriter()
    page = writer.add_blank_page(420, 420)
    page[NameObject('/Resources')] = DictionaryObject()
    stream = DecodedStreamObject()
    stream.set_data(content)
    page[NameObject('/Contents')] = writer._add_object(stream)
    return writer, page


class VectorPathReviewTests(unittest.TestCase):
    def test_blank_page_without_contents_is_a_valid_compression_input(self):
        writer = PdfWriter()
        page = writer.add_blank_page(420, 420)
        compress_pdf_page_streams_losslessly(writer)
        self.assertNotIn('/Contents', page)

    def test_inherited_resources_are_copied_without_losing_or_overwriting_a_form(self):
        writer, page = document(b'/FmM1 Do\n' + (POLYGON + b'\nf\n') * 20)
        existing = DecodedStreamObject()
        existing.set_data(b'1 0 0 rg 5 5 10 10 re f\n')
        reference = writer._add_object(existing)
        resources = DictionaryObject({NameObject('/XObject'): DictionaryObject({
            NameObject('/FmM1'): reference})})
        parent = page['/Parent']
        parent[NameObject('/Resources')] = resources
        del page['/Resources']
        self.assertGreater(deduplicate_pdf_vector_glyphs(writer), 0)
        self.assertIs(parent['/Resources'], resources)
        self.assertEqual(list(resources['/XObject']), ['/FmM1'])
        self.assertEqual(page['/Resources']['/XObject'].raw_get('/FmM1'), reference)
        self.assertIn('/FmM2', page['/Resources']['/XObject'])

    def test_separate_adjusted_strokes_use_outward_inverse_page_bounds(self):
        for line_width in (b'0', b'100'):
            raw = (b'q /Adjust gs ' + line_width + b' w\n'
                   b'2 0 0 -0.5 10 20 cm\n' + POLYGON + b'\nS Q\n') * 20
            writer, page = document(raw)
            page['/Resources'][NameObject('/ExtGState')] = DictionaryObject({
                NameObject('/Adjust'): DictionaryObject({NameObject('/SA'): BooleanObject(True)})})
            self.assertEqual(deduplicate_pdf_vector_glyphs(writer), 1)
            form = next(iter(page['/Resources']['/XObject'].values())).get_object()
            left, bottom, right, top = map(float, form['/BBox'])
            self.assertLessEqual(left, -5)
            self.assertLessEqual(bottom, -800)
            self.assertGreaterEqual(right, 205)
            self.assertGreaterEqual(top, 40)
            self.assertTrue(form.get_data().endswith(b'S\n'))
            self.assertNotIn(b'B\n', form.get_data())

    def test_fill_and_fillstroke_bounds_contain_the_page_at_tiny_scale(self):
        for paint in (b'f', b'f*', b'B', b'B*'):
            raw = (b'q .0001 0 0 .0001 0 0 cm\n' + POLYGON + b'\n' + paint + b' Q\n') * 20
            writer, page = document(raw)
            self.assertEqual(deduplicate_pdf_vector_glyphs(writer), 1)
            form = next(iter(page['/Resources']['/XObject'].values())).get_object()
            left, bottom, right, top = map(float, form['/BBox'])
            self.assertLess(left, 0)
            self.assertLess(bottom, 0)
            self.assertGreater(right, 4200000)
            self.assertGreater(top, 4200000)
            self.assertTrue(form.get_data().endswith(paint + b'\n'))

    def test_large_intermediate_translation_cannot_regain_eligibility_by_cancellation(self):
        prefix = (b'1 0 0 1 10000000000000000 0 cm\n'
                  + b'1 0 0 1 1 0 cm\n' * 64
                  + b'1 0 0 1 -10000000000000000 0 cm\n')
        raw = prefix + (POLYGON + b'\nS\n') * 20
        writer, page = document(raw)
        self.assertEqual(deduplicate_pdf_vector_glyphs(writer), 0)
        self.assertEqual(page.get_contents().get_data(), raw)
        writer, page = document(b'q\n' + raw + b'Q\n' + (POLYGON + b'\nS\n') * 20)
        self.assertEqual(deduplicate_pdf_vector_glyphs(writer), 1)
        data = page.get_contents().get_data()
        self.assertEqual(data.splitlines().count(b'S'), 20)
        self.assertEqual(sum(line.endswith(b' Do') for line in data.splitlines()), 20)

    def test_nonaxis_singular_or_unbounded_stroke_transforms_stay_raw(self):
        for matrix in (b'1 .1 0 1 0 0', b'0 0 0 1 0 0', b'.000000000001 0 0 1 0 0'):
            raw = (b'q ' + matrix + b' cm\n' + POLYGON + b'\nS Q\n') * 20
            writer, page = document(raw)
            self.assertEqual(deduplicate_pdf_vector_glyphs(writer), 0)
            self.assertEqual(page.get_contents().get_data(), raw)

    def test_complete_fill_and_stroke_inherits_state_with_a_shared_safe_bbox(self):
        raw = ((b'q 1 w 10 M\n' + POLYGON + b'\nB Q\n')
               + (b'q 100 w 20 M\n' + POLYGON + b'\nB Q\n')) * 4
        writer, page = document(raw)
        self.assertEqual(deduplicate_pdf_vector_glyphs(writer), 1)
        form = next(iter(page['/Resources']['/XObject'].values())).get_object()
        box = list(map(float, form['/BBox']))
        self.assertLessEqual(box[0], 0)
        self.assertGreaterEqual(box[2], 420)
        self.assertTrue(form.get_data().endswith(b'B\n'))
        self.assertEqual(page.get_contents().get_data().splitlines().count(b'100 w'), 4)

    def test_hairline_adjusted_and_unknown_strokes_stay_in_the_original_stream(self):
        for prefix, state in ((b'0 w ', None),
                              (b'/Adjust gs ', DictionaryObject({NameObject('/SA'): BooleanObject(True)})),
                              (b'/Unknown gs ', None)):
            raw = (prefix + POLYGON + b'\nB\n') * 4
            writer, page = document(raw)
            if state is not None:
                page['/Resources'][NameObject('/ExtGState')] = DictionaryObject({NameObject('/Adjust'): state})
            self.assertEqual(deduplicate_pdf_vector_glyphs(writer), 0)
            self.assertEqual(page.get_contents().get_data(), raw)

    def test_rectangle_bbox_contains_its_true_positive_width_and_height_endpoints(self):
        writer, page = document((RECTANGLE + b'\nf\n') * 20)
        self.assertGreater(deduplicate_pdf_vector_glyphs(writer), 0)
        forms = page['/Resources']['/XObject'].values()
        for reference in forms:
            left, bottom, right, top = map(float, reference.get_object()['/BBox'])
            self.assertLessEqual(left, 80)
            self.assertLessEqual(bottom, 90)
            self.assertGreaterEqual(right, 260)  # 80+180, not the width180
            self.assertGreaterEqual(top, 240)  # 90+150, not the height150

    def test_wide_stroke_is_retained_with_all_inherited_state(self):
        raw = (b'q\n100 w\n10 M\n2 j\n' + POLYGON + b'\nS\nQ\n') * 5
        writer, page = document(raw)
        self.assertEqual(deduplicate_pdf_vector_glyphs(writer), 0)
        self.assertEqual(page.get_contents().get_data(), raw)

    def test_fill_and_stroke_remain_on_their_original_sides_of_state_changes(self):
        raw = (b'q\n/G0 gs\n' + POLYGON + b'\nf\n/G2 gs\n'
               + POLYGON + b'\nS\nQ\n') * 20
        writer, page = document(raw)
        self.assertGreater(deduplicate_pdf_vector_glyphs(writer), 0)
        lines = page.get_contents().get_data().splitlines()
        self.assertEqual(lines.count(b'/G0 gs'), 20)
        self.assertEqual(lines.count(b'/G2 gs'), 20)
        self.assertEqual(lines.count(b'S'), 20)
        self.assertNotIn(b'B', lines)
        for index, line in enumerate(lines):
            if line == b'/G0 gs':
                self.assertTrue(lines[index + 1].endswith(b' Do'))
                self.assertEqual(lines[index + 2], b'/G2 gs')

    def test_no_matching_repetition_leaves_original_bytes_unchanged(self):
        raw = b'q\n' + RECTANGLE + b'\nf\nQ\n'
        writer, page = document(raw)
        self.assertEqual(deduplicate_pdf_vector_glyphs(writer), 0)
        self.assertEqual(page.get_contents().get_data(), raw)

    def test_multiline_literal_path_text_is_not_a_sequence_of_pdf_operators(self):
        raw = (b'BT\n/F1 4 Tf\n5 300 Td\n(path is text\n' + POLYGON
               + b'\nf\n) Tj\nET\n') * 20
        writer, page = document(raw)
        font = DictionaryObject({NameObject('/Type'): NameObject('/Font'),
                                 NameObject('/Subtype'): NameObject('/Type1'),
                                 NameObject('/BaseFont'): NameObject('/Helvetica')})
        page['/Resources'][NameObject('/Font')] = DictionaryObject({
            NameObject('/F1'): writer._add_object(font)})
        self.assertEqual(deduplicate_pdf_vector_glyphs(writer), 0)
        self.assertEqual(page.get_contents().get_data(), raw)

    def test_a_state_change_cannot_separate_an_earlier_pending_subpath_from_its_fill(self):
        triangle = b'10 10 m\n30 10 l\n20 30 l\nh\n'
        raw = (triangle + b'1 0 0 rg\n' + POLYGON + b'\nf\n') * 20
        writer, page = document(raw)
        self.assertEqual(deduplicate_pdf_vector_glyphs(writer), 0)
        self.assertEqual(page.get_contents().get_data(), raw)

    def test_an_existing_form_resource_name_is_not_overwritten(self):
        writer, page = document(b'/FmM1 Do\n' + (POLYGON + b'\nf\n') * 20)
        existing = DecodedStreamObject()
        existing.set_data(b'1 0 0 rg 5 5 10 10 re f\n')
        existing.update({NameObject('/Type'): NameObject('/XObject'),
                         NameObject('/Subtype'): NameObject('/Form'),
                         NameObject('/BBox'): ArrayObject(map(NumberObject, [0, 0, 20, 20]))})
        reference = writer._add_object(existing)
        page['/Resources'][NameObject('/XObject')] = DictionaryObject({NameObject('/FmM1'): reference})
        self.assertGreater(deduplicate_pdf_vector_glyphs(writer), 0)
        objects = page['/Resources']['/XObject']
        self.assertEqual(objects.raw_get('/FmM1'), reference)
        self.assertIn('/FmM2', objects)

    def test_stream_compression_keeps_high_ctm_coordinate_precision_and_text_bytes(self):
        raw = (b'q\n10000 0 0 10000 0 0 cm\n0.010000 0.010000 m\n'
               b'0.019999 0.010000 l\n0.019999 0.020000 l\n'
               b'0.010000 0.020000 l\nh\nf\nQ\n'
               b'BT\n(0.019999 0.010000 m\n) Tj\nET\n')
        writer, page = document(raw)
        compress_pdf_page_streams_losslessly(writer)
        self.assertEqual(page.get_contents().get_data(), raw)

    def test_parsed_operator_serialization_keeps_float_values_and_legal_pdf_number_syntax(self):
        raw = (b'1000000000 0 0 1000000000 -1234567800 -1234567800 cm\n'
               b'1.2345678912345 1.2345678912345 0.00000001 0.00000001 re\nf\n')
        writer, page = document(raw)
        before = ContentStream(page.get_contents(), writer).operations
        serialized = pdf_content_bytes_losslessly(before, writer)
        stream = DecodedStreamObject()
        stream.set_data(serialized)
        after = ContentStream(stream, writer).operations
        self.assertEqual(before, after)
        self.assertIn(b'1.2345678912345', serialized)
        self.assertIn(b'0.00000001', serialized)
        self.assertNotIn(b'e-', serialized)
        self.assertNotIn(b'e+', serialized)


class PDFSpaceAdvanceReviewTests(unittest.TestCase):
    def font(self):
        writer = PdfWriter()
        font = font_resource(writer, space_advance=600)
        return writer, font, font.get_object()['/DescendantFonts'][0].get_object()

    def test_pdf_width_can_override_a_positive_embedded_font_width(self):
        for widths in (ArrayObject([NumberObject(3), ArrayObject([NumberObject(0)])]),
                       ArrayObject([NumberObject(3), NumberObject(3), NumberObject(-1)])):
            with self.subTest(widths=widths):
                writer, font, descendant = self.font()
                descendant[NameObject('/W')] = widths
                self.assertFalse(verified_layout_space_font(font))


    def test_default_and_explicit_pdf_widths_are_resolved_in_order(self):
        writer, font, descendant = self.font()
        descendant[NameObject('/DW')] = NumberObject(0)
        self.assertTrue(verified_layout_space_font(font))  # /W supplies CID3=600
        del descendant['/W']
        self.assertFalse(verified_layout_space_font(font))
        del descendant['/DW']
        self.assertTrue(verified_layout_space_font(font))  # PDF default is1000

    def test_ambiguous_or_malformed_pdf_widths_do_not_prove_a_space(self):
        cases = (
            ArrayObject([NumberObject(3), ArrayObject([NumberObject(600)]),
                         NumberObject(3), NumberObject(3), NumberObject(500)]),
            ArrayObject([NumberObject(3), ArrayObject([BooleanObject(True)])]),
            ArrayObject([NumberObject(3), ArrayObject([FloatObject(float('nan'))])]),
            ArrayObject([NumberObject(3)]),
        )
        for widths in cases:
            with self.subTest(widths=widths):
                writer, font, descendant = self.font()
                descendant[NameObject('/W')] = widths
                self.assertFalse(verified_layout_space_font(font))


class PDFObjectReviewTests(unittest.TestCase):
    def test_decimal_tokens_fit_the_paired_reader_limit_or_fail_before_writing(self):
        for value in (1e-45, -1e-44):
            writer, page = document(b'')
            page[NameObject('/Fixture')] = FloatObject(value)
            preserve_pdf_float_values(writer)
            output = io.BytesIO()
            writer.write(output)
            self.assertEqual(float(PdfReader(io.BytesIO(output.getvalue())).pages[0]['/Fixture']), value)
        for value in (1e-46, -1e-45, 1e-210, 1e200):
            writer, page = document(b'')
            page[NameObject('/Fixture')] = FloatObject(value)
            preserve_pdf_float_values(writer)
            with self.assertRaisesRegex(ValueError, 'numeric-token limit'):
                writer.write(io.BytesIO())

    def test_container_precision_keeps_stream_and_indirect_identities(self):
        writer, page = document(b'0 0 20 20 re f')
        reference = writer._add_object(FloatObject(1.2345678912345))
        payload = page.raw_get('/Contents').get_object()
        array = ArrayObject([reference, FloatObject(0.00000000012345678912345)])
        page[NameObject('/Fixture')] = DictionaryObject({NameObject('/Values'): array})
        preserve_pdf_float_values(writer)
        self.assertIs(page['/Fixture']['/Values'], array)
        self.assertIs(array[0], reference)
        self.assertIs(page.raw_get('/Contents').get_object(), payload)
        for custom in (False, True):
            output = io.BytesIO()
            if custom:
                write_pdf_with_object_streams(writer, output)
            else:
                writer.write(output)
            observed = PdfReader(io.BytesIO(output.getvalue())).pages[0]['/Fixture']['/Values']
            self.assertEqual(float(observed[0].get_object()), 1.2345678912345)
            self.assertEqual(float(observed[1]), 0.00000000012345678912345)

    def test_marked_is_a_parsed_pdf_boolean_value(self):
        for value, expected in ((BooleanObject(True), True), (BooleanObject(False), False),
                                (NumberObject(7), False), (NameObject('/True'), False)):
            writer, _ = document(b'')
            writer.root_object[NameObject('/MarkInfo')] = DictionaryObject({NameObject('/Marked'): value})
            output = io.BytesIO()
            writer.write(output)
            reader = PdfReader(io.BytesIO(output.getvalue()))
            self.assertIs(pdf_marked_true(reader), expected)

    def test_extgstate_font_does_not_leave_a_visible_glyph_marked_as_decoration(self):
        raw = b'BT /Fblank 12 Tf q /Visible gs <0003> Tj Q <0003> Tj ET'
        writer, page, _, _ = page_fixture(raw)
        visible = font_resource(writer, visible_space=True)
        page['/Resources'][NameObject('/ExtGState')] = DictionaryObject({
            NameObject('/Visible'): DictionaryObject({
                NameObject('/Font'): ArrayObject([visible, NumberObject(12)])})})
        tag_pdf_formulas(writer, [])
        self.assertEqual(artifact_texts(ContentStream(page.get_contents(), writer).operations), [b'\x00\x03'])



if __name__ == '__main__':
    unittest.main()
