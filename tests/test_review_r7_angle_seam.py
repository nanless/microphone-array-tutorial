"""R363 principal-seam angle differences and their public consumers."""
import math
import unittest
import numpy as np

from codes.chapters.ch09.core.tracking import angle_difference, ConstantVelocityKalman, CircularParticleFilter
from codes.chapters.ch09.core.echo_av_fusion import (
    small_set_distances, angular_ospa_distance, clear_mot_angular_summary,
    echo_sector_gate, fov_bernoulli_missed_detection_update, oosm_cv_retrodiction_update,
)


class PrincipalSeamTests(unittest.TestCase):
    def test_shared_difference_and_two_ospa_entries(self):
        first = float(np.nextafter(180., 0.))
        gap = 2.**-45  # Exact circular distance from first to -180 degrees.
        for a, b, expected in ((-180., first, gap), (first, -180., -gap),
                               (100., np.nextafter(100., math.inf), -2.**-46),
                               (0., 1e-200, -1e-200), (180., -180., 0.)):
            self.assertEqual(angle_difference(a, b), expected)
            for order in (1, 2):
                self.assertEqual(angular_ospa_distance([a], [b], order=order), abs(expected))
                self.assertEqual(small_set_distances([a], [b], order=order)['ospa'], abs(expected))
        for a, b in ((0., 180.), (180., 0.), (90., -90.), (-90., 90.)):
            self.assertEqual(angle_difference(a, b), -180.)
        for order in (1, 2):
            self.assertEqual(small_set_distances([first, -180.], [-180., first], order=order)['ospa'], 0.)
        self.assertEqual(angle_difference(np.finfo(float).max, np.finfo(float).max), 0.)

    def test_signed_innovation_and_particle_posterior(self):
        first = float(np.nextafter(180., 0.)); gap = 2.**-45
        result = oosm_cv_retrodiction_update([first, 0.], np.eye(2), -180., 1., 0.)
        self.assertEqual(result['retrodiction']['innovation_deg'], gap)
        tracker = ConstantVelocityKalman([first, 0.], np.eye(2), np.zeros((2, 2)))
        self.assertEqual(tracker.update(-180., 1.)[0], -180.)
        particles = CircularParticleFilter([first, -180.])
        particles.update(-180., gap/2, clutter_probability=0.)
        np.testing.assert_allclose(particles.weights,
                                   [1/(1+math.exp(2)), 1/(1+math.exp(-2))], rtol=1e-14, atol=0.)

    def test_small_legal_gates_and_clear_matching(self):
        first = float(np.nextafter(180., 0.)); gap = 2.**-45
        sector = echo_sector_gate([first], speaker_angle_deg=-180.,
                                  sector_half_width_deg=gap/2, far_end_active=True, double_talk=False)
        self.assertEqual(sector['angular_distance_to_speaker_deg'][0], gap)
        self.assertFalse(sector['inside_sector_mask'][0])
        fov = fov_bernoulli_missed_detection_update(.8, first, camera_center_deg=-180.,
                                                    fov_half_width_deg=gap/2, in_fov_detection_prob=.9)
        self.assertFalse(fov['in_fov'])
        self.assertEqual(fov['fov_aware_posterior_existence'], .8)
        result = clear_mot_angular_summary([{1:first, 2:-180.}], [{1:-180., 2:first}])
        self.assertEqual(result['frames'][0]['matches'], {'1':2, '2':1})
        self.assertEqual(result['motp_deg'], 0.)
        carry = clear_mot_angular_summary([{1:first}, {1:first}], [{7:first}, {7:-180.}],
                                          match_threshold_deg=gap/2)
        self.assertEqual(carry['frames'][1]['matches'], {})
        self.assertEqual((carry['fn'], carry['fp']), (1, 1))

    def test_near_antipodes_keep_direction_and_slew_step(self):
        from codes.chapters.ch09.tracking_crossing_dropout_demo import _limited_step
        near_half_turn = float(np.nextafter(180., 0.))
        for epsilon in (1e-14, 1e-200):
            self.assertEqual(angle_difference(-epsilon, -180.), near_half_turn)
            self.assertEqual(angle_difference(-180., -epsilon), -near_half_turn)
            self.assertEqual(angle_difference(epsilon, -180.), -near_half_turn)
            self.assertEqual(angle_difference(-180., epsilon), near_half_turn)
            self.assertEqual(_limited_step(-180., -epsilon, 5., 1.), -175.)
            self.assertEqual(_limited_step(-180., epsilon, 5., 1.), 175.)
        # A rotated near-antipode and genuine half-turn controls.
        other = float(np.nextafter(-90., 0.))
        self.assertEqual(angle_difference(90., other), near_half_turn)
        self.assertEqual(angle_difference(other, 90.), -near_half_turn)
        self.assertEqual(angle_difference(90., -90.), -180.)
        self.assertEqual(angle_difference(-90., 90.), -180.)
