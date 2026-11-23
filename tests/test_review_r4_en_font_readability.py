"""Independent font/overflow controls for the reviewed English figure labels."""
import unittest
from unittest.mock import patch

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from scripts import make_figures_en as en


class EnglishFigureReadabilityTests(unittest.TestCase):
    def tearDown(self):
        plt.close('all')

    def flowchart(self):
        fig, axes = plt.subplots(2, 2, figsize=(9.5, 9))
        ax = fig.axes[0]
        ax.set(xlim=(0, 10), ylim=(0, 10))
        ax.axis('off')
        for y, text in ((8, 'Playback\n$x[n]$'),
                        (3, 'NLMS\n$\\hat{\\vec{w}}[n+1]$')):
            ax.add_patch(Rectangle((2, y - 1), 4, 2))
            ax.text(4, y, text, ha='center', va='center')
        note_axis = fig.axes[3]
        note_axis.set(xlim=(0, 1), ylim=(0, 1))
        note_axis.axvspan(.49, .51, alpha=.1)
        note = note_axis.text(.5, .5, 'Synthetic\nnear-end burst',
                              ha='center', va='center')
        protected = en._publication_font_artists(fig, 'fig18_aec.png')
        return fig, note, protected

    def test_highlight_is_not_an_annotation_container(self):
        fig, note, protected = self.flowchart()
        en._autofit_translated_layout(fig, {id(note)}, protected)
        en._validate_publication_fonts(fig, 'fig18_aec.png', protected)
        self.assertEqual(note.get_fontsize(), 11)
        span = fig.axes[3].patches[0].get_window_extent(fig.canvas.get_renderer())
        self.assertGreater(note.get_window_extent(fig.canvas.get_renderer()).width,
                           span.width)

    def test_small_font_and_outside_page_annotation_are_rejected(self):
        for failure in ('small', 'outside'):
            with self.subTest(failure=failure):
                fig, note, protected = self.flowchart()
                if failure == 'small':
                    note.set_fontsize(10)
                else:
                    note.set_position((10, .5))
                with self.assertRaisesRegex(ValueError, 'annotation'):
                    en._validate_publication_fonts(fig, 'fig18_aec.png', protected)

    def test_flowchart_box_overflow_is_rejected(self):
        fig, _, protected = self.flowchart()
        # Keep the original protected artist identity, but exceed its actual box.
        label = fig.axes[0].texts[0]
        label.set_text('Playback reference with excessive width')
        with self.assertRaisesRegex(ValueError, 'flowchart label'):
            en._validate_publication_fonts(fig, 'fig18_aec.png', protected)

    def test_full_size_title_and_legend_have_independent_overflow_controls(self):
        fig, axes = plt.subplots(1, 3, figsize=(9.5, 4.3))
        for ax in axes:
            ax.set_title('Short\ntitle')
            ax.plot([0, 1], [0, 1], label='Short')
            ax.legend()
        protected = en._publication_font_artists(fig, 'fig61_equal_drr_response.png')
        en._autofit_translated_layout(fig, set(), protected)
        en._validate_publication_fonts(fig, 'fig61_equal_drr_response.png', protected)
        for artist in [axes[0].title, axes[2].get_legend().get_texts()[0]]:
            with self.subTest(role=type(artist).__name__):
                original = artist.get_text()
                artist.set_text('Excessive width ' * 20)
                with self.assertRaisesRegex(ValueError, 'full-size'):
                    en._validate_publication_fonts(fig, 'fig61_equal_drr_response.png', protected)
                artist.set_text(original)
        axes[2].get_legend().get_texts()[0].set_fontsize(10)
        with self.assertRaisesRegex(ValueError, 'full-size legend'):
            en._validate_publication_fonts(fig, 'fig61_equal_drr_response.png', protected)

    def test_figure_local_translation_does_not_change_other_consumers(self):
        translator = en.FigureTranslator(en.TRANSLATIONS_PATH)
        fig, ax = plt.subplots()
        text = ax.text(.5, .5, 'High-tone inverted')
        en._apply_figure_text_overrides(fig, 'unrelated.png', translator.figure_overrides)
        self.assertEqual(text.get_text(), 'High-tone inverted')
        with patch.object(en, 'PRIMARY_CHARMAP', set(range(128))):
            en._apply_figure_text_overrides(fig, 'probe.png',
                                           {'probe.png': {'High-tone inverted': 'High tone: 180°'}})
        self.assertEqual(text.get_text(), 'High tone: 180°')
        self.assertEqual(text.get_fontfamily(), ['DejaVu Sans'])
        with self.assertRaisesRegex(ValueError, 'missing'):
            en._apply_figure_text_overrides(fig, 'probe.png',
                                           {'probe.png': {'Unseen label': 'Short'}})

if __name__ == '__main__':
    unittest.main()
