"""Independent fraction and analytic expectations for E06-08..E06-10."""
import json
from pathlib import Path
import subprocess
import sys
import unittest

from codes.chapters.ch06.aec_advanced_exercises import run_exercises


class AdvancedAECExercisesTest(unittest.TestCase):
    def test_stable_ids_and_json(self):
        result = run_exercises()
        self.assertEqual(list(result), ["E06-08", "E06-09", "E06-10"])
        json.dumps(result, allow_nan=False)
        script = Path(__file__).resolve().parents[1] / "codes/chapters/ch06/aec_advanced_exercises.py"
        proc = subprocess.run([sys.executable, str(script)], capture_output=True, text=True, check=True)
        self.assertEqual(json.loads(proc.stdout), result)

    def test_ipnlms_zero_start_deadlock(self):
        self.assertEqual(run_exercises()["E06-08"]["pure_proportionate_zero_start_update"], [0, 0])

    def test_rls_forgetting_time_conversion(self):
        result = run_exercises()
        self.assertAlmostEqual(result["E06-09"]["old_observation_weight_after_100_updates"], .3660323412732292)
        self.assertEqual(result["E06-09"]["elapsed_seconds_if_16khz_sample_updates"], .00625)
        self.assertEqual(result["E06-09"]["elapsed_seconds_if_10ms_block_updates"], 1.)

    def test_two_band_haar_block_cost(self):
        result = run_exercises()
        self.assertEqual(result["E06-10"]["physical_taps"], 4)
        self.assertEqual(result["E06-10"]["minimum_block_taps"], 3)
        self.assertEqual(result["E06-10"]["unconstrained_crossband_weights"], 12)
        self.assertEqual(result["E06-10"]["minimum_pairing_wait_samples"], 1)


if __name__ == "__main__":
    unittest.main()
