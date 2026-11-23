"""Independent counterexamples retained from the multiround chapter review."""
import math
import copy
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from codes.chapters.ch10.core.pdm_and_eval import (
    cic_decimation_metrics, p340_and_p808_evaluation,
)
from codes.chapters.ch12.core.functional_beamforming import functional_beamforming_map
from codes.chapters.ch13.core.dasf_and_ptp import centralized_lcmv
from codes.chapters.ch14.core.ambisonic_sofa_rendering import audit_sofa_simple_freefield_geometry
from codes.chapters.ch14.examples import generate_wearable_audio
from codes.chapters.ch15.core.av_alignment import map_video_to_audio_frames
from codes.chapters.ch15.examples import generate_tse_audio


class EngineeringReviewTests(unittest.TestCase):
    def test_zero_attenuation_does_not_establish_joint_duplex_class(self):
        result = p340_and_p808_evaluation(
            (0.2, 0.2), {"bypass": (0.2, 0.2)},
            [{"worker_id": "w", "trap": 1, "gold": (5, 1), "clips": [4]}],
        )
        row = result["p340"]["bypass"]
        self.assertEqual(row["attenuation_type"], "Type 1")
        self.assertIsNone(row["joint_type"])
        self.assertEqual(row["joint_status"], "insufficient_evidence")
        self.assertIn("double_talk_telr", row["missing_evidence"])

    def test_cic_comb_delay_changes_filter_but_not_output_nyquist(self):
        f_clock, r, m, order, frequency = 3_072_000, 16, 2, 3, 60_000
        result = cic_decimation_metrics(f_clock, r, order, frequency, diff_delay=m)
        # Construct the equivalent FIR directly, independently of the sine ratio.
        impulse = np.ones(r * m) / (r * m)
        fir = np.array([1.0])
        for _ in range(order):
            fir = np.convolve(fir, impulse)
        response = abs(np.sum(fir * np.exp(-2j * math.pi * frequency / f_clock * np.arange(len(fir)))))
        self.assertEqual(result["output_rate_hz"], 192_000)
        self.assertAlmostEqual(result["passband_gain"], response, places=14)
        with self.assertRaises(ValueError):
            cic_decimation_metrics(f_clock, r, order, 96_000, diff_delay=m)


class ImagingDistributedReviewTests(unittest.TestCase):
    def test_nu_one_retains_small_positive_power(self):
        result = functional_beamforming_map(np.diag([1.0, 1e-13]), np.eye(2), nu=1)
        np.testing.assert_allclose(result, [1.0, 1e-13], rtol=1e-14, atol=0)

    def test_lcmv_independent_column_units_leave_solution_unchanged(self):
        # Square independent constraints determine w directly, regardless of R.
        cov = np.array([[2, 0.3j], [-0.3j, 1]])
        base = np.array([[1, 1j], [1j, 2]], dtype=complex)
        expected = np.array([1 + 0.2j, -0.4j])
        response = base.conj().T @ expected
        for scales in ([1, 1e-8], [1e-100, 1e100]):
            c = base * scales
            f = response * scales
            actual = centralized_lcmv(cov, c, f)["weights"]
            np.testing.assert_allclose(actual, expected, rtol=2e-14, atol=2e-14)
        with self.assertRaises(ValueError):
            centralized_lcmv(np.eye(2), np.ones((2, 2)), [1, 1])


    def test_lcmv_raw_diagnostic_does_not_overflow_before_exact_cancellation(self):
        constraints = 1e308 * (np.ones((3, 3)) - 2*np.eye(3))
        result = centralized_lcmv(np.eye(3), constraints, np.full(3, 1e308))
        np.testing.assert_array_equal(result["weights"], np.ones(3))
        self.assertEqual(result["constraint_max_abs_residual"], 0.)
        self.assertEqual(result["noise_power"], 3.)

    def test_lcmv_rejects_response_lost_during_column_normalization(self):
        for scale, response in ((1e308, 1e-20), (1e200, 1e-200)):
            for value in (response, 1j*response, 1 + 1j*response):
                with self.subTest(scale=scale, value=value):
                    with self.assertRaisesRegex(ValueError, "normalized response is below"):
                        centralized_lcmv(np.eye(2), scale*np.eye(2), np.full(2, value))
        # A truly zero response stays feasible and must not be mistaken for loss.
        result = centralized_lcmv(np.eye(2), 1e308*np.eye(2), np.zeros(2))
        np.testing.assert_array_equal(result["weights"], np.zeros(2))
        self.assertEqual(result["constraint_max_abs_residual"], 0.)

    def test_lcmv_diagnostic_rejects_unrepresentable_original_units(self):
        # Isolate report validation from the already-checked unique solver.
        module = "codes.chapters.ch13.core.dasf_and_ptp.lcmv_weights"
        constraints = 1e308 * (np.ones((3, 3)) - 2*np.eye(3))
        with patch(module, return_value=np.full(3, 4.)):
            with self.assertRaisesRegex(ValueError, "diagnostics exceed"):
                centralized_lcmv(np.eye(3), constraints, np.full(3, 1e308))
        with patch(module, return_value=np.full(2, 1e-200)):
            with self.assertRaisesRegex(ValueError, "below float64 support"):
                centralized_lcmv(np.eye(2), 1e-200*np.eye(2), np.zeros(2))
        with self.assertRaisesRegex(ValueError, "positive LCMV noise power"):
            centralized_lcmv(np.eye(2), np.eye(2), np.full(2, 1e-200))


class WearableTseReviewTests(unittest.TestCase):
    def test_listener_yaw_does_not_rotate_local_receiver_coordinates_twice(self):
        ears = np.array([[0, 0.09, 0], [0, -0.09, 0]])
        result = audit_sofa_simple_freefield_geometry([0, 1, 0], [0, 0, 1], ears)
        self.assertTrue(result["passed"])
        self.assertFalse(audit_sofa_simple_freefield_geometry([0, 1, 0], [0, 0, 1], ears[::-1])["ear_order_valid"])

    def test_online_mapping_matches_independent_arrival_enumeration(self):
        for delay in (0, 0.04, 0.08, 0.12):
            result = map_video_to_audio_frames(30, 0.016, 25, delay)
            arrivals = (np.arange(100) + 0.5) / 25 + delay
            expected = np.searchsorted(arrivals, result["audio_center_times_s"], side="right") - 1
            np.testing.assert_array_equal(result["video_frame_indices"], expected)
            valid = result["valid_video_mask"]
            self.assertTrue(np.all(result["mapped_video_arrival_times_s"][valid] <= result["audio_center_times_s"][valid]))
            self.assertEqual(result["video_frame_indices"][0], -1)
        boundary = map_video_to_audio_frames(1, 0.04, 25)
        self.assertTrue(boundary["valid_video_mask"][0])
        self.assertEqual(boundary["mapped_video_arrival_times_s"][0], 0.02)
        offline = map_video_to_audio_frames(15, 0.016, 25, alignment_mode="offline_nearest_exposure")
        self.assertEqual(offline["video_frame_indices"][0], 0)
        self.assertGreater(offline["mapped_video_arrival_times_s"][0], offline["audio_center_times_s"][0])

    def test_analytic_corruption_is_rejected_before_first_directory_creation(self):
        for generator, key, field in (
            (generate_wearable_audio, "wearable_fused_output", "residual_noise_power"),
            (generate_wearable_audio, "wearable_fused_output", "reference_mse_ch0"),
            (generate_tse_audio, "bubble_d10cm", "sensor_noise_output_power"),
            (generate_tse_audio, "bubble_d10cm", "reference_mse"),
            (generate_tse_audio, "bubble_d10cm", "reference_nmse"),
        ):
            report, arrays = generator.run_experiment()
            for mutation in (99.0, True, None, float("nan"), -1e-16):
                corrupt = copy.deepcopy(report)
                if mutation is None:
                    del corrupt["analytic_expectations"][key][field]
                else:
                    corrupt["analytic_expectations"][key][field] = mutation
                with tempfile.TemporaryDirectory() as td:
                    target = Path(td) / "not-created"
                    with patch.object(generator, "run_experiment", return_value=(corrupt, arrays)), patch.object(Path, "mkdir") as mkdir:
                        with self.assertRaisesRegex(ValueError, "independent analytic"):
                            generator.generate_assets(target)
                        mkdir.assert_not_called()
                    self.assertFalse(target.exists())


if __name__ == "__main__":
    unittest.main()
