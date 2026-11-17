"""Independent actual-PCM publication and strict numeric-type negative controls for Chapter 16."""
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
WEARABLE_SPECS = (
    (81, "fig81_foa_intensity_diffuseness.png", "figure81_foa_intensity_diffuseness.json"),
    (82, "fig82_bsm_ls_vs_magls.png", "figure82_bsm_ls_vs_magls.json"),
    (83, "fig83_own_voice_fusion.png", "figure83_own_voice_fusion.json"),
    (84, "fig84_hear_through_comb.png", "figure84_hear_through_comb.json"),
)


class WearablePublicationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.audio = self.root / "codes/chapters/ch16/wearable_audio"
        shutil.copytree(ROOT / "codes/chapters/ch16/wearable_audio", self.audio)
        self.manifest = json.loads((self.audio / "MANIFEST.json").read_bytes())

        report_dir = self.root / "codes/chapters/ch16/reports"
        report_dir.mkdir(parents=True)
        (self.root / "figures").mkdir()
        for _num, png_name, json_name in WEARABLE_SPECS:
            shutil.copy2(ROOT / "codes/chapters/ch16/reports" / json_name, report_dir / json_name)
            shutil.copy2(ROOT / "figures" / png_name, self.root / "figures" / png_name)

        report83 = json.loads((report_dir / "figure83_own_voice_fusion.json").read_bytes())
        for relative in report83["source_sha256"]:
            target = self.root / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / relative, target)

        self.site = self.root / "site"
        shutil.copytree(self.audio, self.site / "wearable_audio")
        (self.site / "research").mkdir()
        wav_names = sorted(qc.WEARABLE_AUDIO_CHANNELS)
        for page, prefix in (
            (self.site / "16_spatial-wearable-arrays.html", ""),
            (self.site / "research/05_exercises_and_audio.html", "../"),
        ):
            page.write_text(
                "\n".join(
                    f'<audio src="{prefix}wearable_audio/{name}" controls preload="none" aria-label="wearable"></audio>'
                    for name in wav_names
                )
                + f'<a href="{prefix}wearable_audio/MANIFEST.json">Manifest</a>',
                encoding="utf-8",
            )

    def check(self, function):
        errors = []
        with (
            patch.object(qc, "ROOT", self.root),
            patch.object(qc, "WEARABLE_AUDIO_ROOT", self.audio),
            patch.object(qc, "SITE", self.site),
        ):
            function(errors)
        return errors

    def test_actual_audio_copies_and_independent_integer_scores_pass(self):
        self.assertEqual(self.check(qc.check_wearable_audio), [])
        self.assertEqual(self.check(qc.check_wearable_figures), [])

    def test_extra_file_and_changed_pcm_fail_without_rewriting(self):
        target = self.site / "wearable_audio/wearable_fused_output.wav"
        original = target.read_bytes()
        target.write_bytes(original[:-2] + (b"\x01\x00" if original[-2:] != b"\x01\x00" else b"\x02\x00"))
        changed = target.read_bytes()
        self.assertTrue(self.check(qc.check_wearable_audio))
        self.assertEqual(target.read_bytes(), changed)
        target.write_bytes(original)
        (self.site / "wearable_audio/unlisted").write_bytes(b"not an asset")
        self.assertTrue(self.check(qc.check_wearable_audio))

    def test_figure_rejects_equal_valued_wrong_types_with_matching_png_report_sha(self):
        report_path = self.root / "codes/chapters/ch16/reports/figure83_own_voice_fusion.json"
        original = json.loads(report_path.read_bytes())
        real_open = Image.open
        with real_open(self.root / "figures/fig83_own_voice_fusion.png") as picture:
            info83 = dict(picture.info)

        for field in ("figure", "mse_bool", "bad_exercise"):
            with self.subTest(field=field):
                report = copy.deepcopy(original)
                if field == "figure":
                    report["figure"] = 83.0
                elif field == "mse_bool":
                    report["pcm_reference_mse_ch0"]["wearable_fused_output"] = False
                else:
                    report["exercise_ids"] = [16]
                report_path.write_text(json.dumps(report, allow_nan=False), encoding="utf-8")
                updated_digest = hashlib.sha256(report_path.read_bytes()).hexdigest()

                def fake_open(path, *args, **kwargs):
                    if Path(path).name == "fig83_own_voice_fusion.png":
                        metadata = {**info83, "NumericalReportDigest": updated_digest}

                        class Picture:
                            def __enter__(self_inner):
                                self_inner.info = metadata
                                return self_inner

                            def __exit__(self_inner, *exc):
                                return False

                        return Picture()
                    return real_open(path, *args, **kwargs)

                with patch("PIL.Image.open", side_effect=fake_open):
                    errors = self.check(qc.check_wearable_figures)
                self.assertTrue(errors)
                self.assertIn("actual audio/report mismatch", "\n".join(errors))

    def test_missing_and_hidden_players_fail(self):
        path = self.site / "16_spatial-wearable-arrays.html"
        text = path.read_text(encoding="utf-8")
        path.write_text(text.replace("<audio", "<audio hidden", 1), encoding="utf-8")
        self.assertTrue(self.check(qc.check_wearable_audio))
        path.write_text(text.replace("controls", "autoplay", 1), encoding="utf-8")
        self.assertTrue(self.check(qc.check_wearable_audio))


if __name__ == "__main__":
    unittest.main()
