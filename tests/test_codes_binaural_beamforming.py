import unittest
import wave

import numpy as np

from codes.chapters.ch05.examples import binaural_beamforming_demo as demo


class BinauralBeamformingTest(unittest.TestCase):
    def test_matched_interference_power_does_not_imply_matched_spatial_cues(self):
        result = demo.run_exercises()
        directional = result["E05-25"]
        np.testing.assert_allclose(directional["bmvdr"], [1/3, 1/3])
        np.testing.assert_allclose(directional["blcmv"], [1/3, -1/3])
        np.testing.assert_allclose(directional["null"], [0, 0], atol=1e-15)
        diffuse = result["E05-26"]["blcmv"]
        np.testing.assert_allclose(diffuse, [[1/3, 1/3], [1/3, 1/2]])
        self.assertAlmostEqual(diffuse[0][1] / np.sqrt(diffuse[0][0] * diffuse[1][1]),
                               np.sqrt(2/3))
        np.testing.assert_allclose(
            result["E05-27"]["complex_right_response_real_imag"], [0, 1/3], atol=1e-15)

    def test_actual_stereo_pcm_preserves_distinct_signed_tone_and_integer_scores(self):
        manifest = demo.check_assets()
        n = np.arange(2400, 29600)
        expected = {"input": [.12, -.12], "bmvdr": [.04, .04],
                    "blcmv": [.04, -.04], "null": [0., 0.]}
        for role, interferer in expected.items():
            row = manifest["measurements"][role]
            with wave.open(str(demo.OUTPUT / f"binaural_{role}.wav"), "rb") as source:
                self.assertEqual((source.getframerate(), source.getnchannels(),
                                  source.getsampwidth(), source.getnframes()),
                                 (16000, 2, 2, 32000))
                q = np.frombuffer(source.readframes(32000), "<i2").reshape(-1, 2)[n].astype(np.int64)
            phasors = [
                (2 / len(n)) * (np.exp(-2j * np.pi * hz * n / 16000) @ q) / 32768
                for hz in [500, 1500]
            ]
            np.testing.assert_allclose(phasors, [[.12, .12], interferer], atol=2e-5, rtol=0)
            score = row["actual_pcm_total"]
            self.assertEqual(score["integer_energy_each_ear"], np.sum(q**2, axis=0).tolist())
            self.assertEqual(score["integer_cross_sum"], int(np.sum(q[:, 0] * q[:, 1])))
            self.assertEqual(score["integer_mean_square_denominator_per_ear"], 27200 * 32768**2)
