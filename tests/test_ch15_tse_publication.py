"""Independent actual-PCM publication and strict numeric-type negative controls for Chapter 15."""
import copy
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from PIL import Image
from scripts import quality_check as qc

ROOT = Path(__file__).resolve().parents[1]
TSE_SPECS = (
    (80, "fig80_tse_conditioning_modes.png", "figure80_tse_conditioning_modes.json"),
    (81, "fig81_nearfield_distance_bubble.png", "figure81_nearfield_distance_bubble.json"),
    (82, "fig82_diarization_tsvad_gss.png", "figure82_diarization_tsvad_gss.json"),
    (83, "fig83_av_asynchrony_occlusion.png", "figure83_av_asynchrony_occlusion.json"),
)


class TsePublicationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.audio = self.root / "codes/chapters/ch15/tse_audio"
        shutil.copytree(ROOT / "codes/chapters/ch15/tse_audio", self.audio)
        self.manifest = json.loads((self.audio / "MANIFEST.json").read_bytes())

        report_dir = self.root / "codes/chapters/ch15/reports"
        report_dir.mkdir(parents=True)
        (self.root / "figures").mkdir()
        for _num, png_name, json_name in TSE_SPECS:
            shutil.copy2(ROOT / "codes/chapters/ch15/reports" / json_name, report_dir / json_name)
            shutil.copy2(ROOT / "figures" / png_name, self.root / "figures" / png_name)

        report81 = json.loads((report_dir / "figure81_nearfield_distance_bubble.json").read_bytes())
        for relative in report81["source_sha256"]:
            target = self.root / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / relative, target)

        self.site = self.root / "site"
        shutil.copytree(self.audio, self.site / "tse_audio")
        (self.site / "research").mkdir()
        wav_names = sorted(qc.TSE_AUDIO_CHANNELS)
        for page, prefix in (
            (self.site / "15_target-extraction-diarization.html", ""),
            (self.site / "research/05_exercises_and_audio.html", "../"),
        ):
            page.write_text(
                "\n".join(
                    f'<audio src="{prefix}tse_audio/{name}" controls preload="none" aria-label="tse"></audio>'
                    for name in wav_names
                )
                + f'<a href="{prefix}tse_audio/MANIFEST.json">Manifest</a>',
                encoding="utf-8",
            )

    def check(self, function):
        errors = []
        with (
            patch.object(qc, "ROOT", self.root),
            patch.object(qc, "TSE_AUDIO_ROOT", self.audio),
            patch.object(qc, "SITE", self.site),
        ):
            function(errors)
        return errors

    def test_actual_audio_copies_and_independent_integer_scores_pass(self):
        self.assertEqual(self.check(qc.check_tse_audio), [])
        self.assertEqual(self.check(qc.check_tse_figures), [])

    def test_extra_file_and_changed_pcm_fail_without_rewriting(self):
        target = self.site / "tse_audio/tse_bubble_d10cm.wav"
        original = target.read_bytes()
        target.write_bytes(original[:-2] + (b"\x01\x00" if original[-2:] != b"\x01\x00" else b"\x02\x00"))
        changed = target.read_bytes()
        self.assertTrue(self.check(qc.check_tse_audio))
        self.assertEqual(target.read_bytes(), changed)
        target.write_bytes(original)
        (self.site / "tse_audio/unlisted").write_bytes(b"not an asset")
        self.assertTrue(self.check(qc.check_tse_audio))

    def test_figure_rejects_equal_valued_wrong_types_with_matching_png_report_sha(self):
        report_path = self.root / "codes/chapters/ch15/reports/figure81_nearfield_distance_bubble.json"
        original = json.loads(report_path.read_bytes())
        real_open = Image.open
        with real_open(self.root / "figures/fig81_nearfield_distance_bubble.png") as picture:
            info81 = dict(picture.info)

        for field in ("figure", "nmse_bool", "bad_exercise"):
            with self.subTest(field=field):
                report = copy.deepcopy(original)
                if field == "figure":
                    report["figure"] = 81.0
                elif field == "nmse_bool":
                    report["pcm_integer_nmse"]["bubble_d10cm"] = False
                else:
                    report["exercise_ids"] = [17]
                report_path.write_text(json.dumps(report, allow_nan=False), encoding="utf-8")
                updated_digest = hashlib.sha256(report_path.read_bytes()).hexdigest()

                def fake_open(path, *args, **kwargs):
                    if Path(path).name == "fig81_nearfield_distance_bubble.png":
                        metadata = {**info81, "NumericalReportDigest": updated_digest}

                        class Picture:
                            def __enter__(self_inner):
                                self_inner.info = metadata
                                return self_inner

                            def __exit__(self_inner, *exc):
                                return False

                        return Picture()
                    return real_open(path, *args, **kwargs)

                with patch("PIL.Image.open", side_effect=fake_open):
                    errors = self.check(qc.check_tse_figures)
                self.assertTrue(errors)
                self.assertIn("actual audio/report mismatch", "\n".join(errors))

    def test_missing_and_hidden_players_fail(self):
        path = self.site / "15_target-extraction-diarization.html"
        text = path.read_text(encoding="utf-8")
        path.write_text(text.replace("<audio", "<audio hidden", 1), encoding="utf-8")
        self.assertTrue(self.check(qc.check_tse_audio))
        path.write_text(text.replace("controls", "autoplay", 1), encoding="utf-8")
        self.assertTrue(self.check(qc.check_tse_audio))


if __name__ == "__main__":
    unittest.main()
