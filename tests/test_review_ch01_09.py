"""Independent numerical regressions from the 2026-11-23 chapter review."""

from decimal import Decimal, localcontext
import math
import unittest
import warnings

import numpy as np

from codes.chapters.ch02.core.diffuse_and_fractional import bessel_j0
from codes.chapters.ch03.core.baffle_and_vent import acoustic_vent_highpass
from codes.chapters.ch08.core.separation_metrics_streaming import (
    causal_fir_toeplitz, ci_sdr_toeplitz, single_speaker_bypass_gate,
)


def decimal_j0(value):
    """Independent high-precision series; never calls the production helper."""
    with localcontext() as context:
        context.prec = 100
        quarter_square = Decimal.from_float(float(value)) ** 2 / 4
        term = Decimal(1)
        total = term
        for index in range(1, 600):
            term *= -quarter_square / (index * index)
            total += term
            if abs(term) < Decimal("1e-90"):
                return float(total)
    raise AssertionError("independent J0 series did not converge")


class SpecialFunctionReviewTests(unittest.TestCase):
    def test_j0_cancellation_and_allowed_series_budget(self):
        expected = 0.047689310796833535
        self.assertEqual(decimal_j0(12.0), expected)
        for budget in (20, 80):
            self.assertAlmostEqual(float(bessel_j0(12.0, max_terms=budget)),
                                   expected, delta=2e-15)

    def test_j0_branch_boundaries_against_decimal(self):
        values = np.array([0.0, 1.0, 4.0, np.nextafter(4.0, math.inf),
                           12.0, 56.0, np.nextafter(56.0, math.inf), 64.0, 100.0])
        expected = np.array([decimal_j0(value) for value in values])
        np.testing.assert_allclose(bessel_j0(values), expected, rtol=0, atol=3e-15)
        np.testing.assert_array_equal(bessel_j0(-values), bessel_j0(values))

    def test_j0_large_finite_arguments_do_not_overflow(self):
        values = np.array([1e100, 1e308, np.finfo(float).max])
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            actual = bessel_j0(values)
        self.assertTrue(np.all(np.isfinite(actual)))
        # The leading Hankel amplitude bounds these extremely large values;
        # this checks range handling, not a universal phase-accuracy claim.
        amplitude = math.sqrt(2 / math.pi) / np.sqrt(values)
        self.assertTrue(np.all(np.abs(actual) <= amplitude * (1 + 1e-14)))
        self.assertTrue(np.all(actual != 0))

    def test_j0_near_hankel_boundary_keeps_approximation_error_explicit(self):
        # Independent review found a larger omitted-term error here than at
        # the other selected points. This is a finite regression, not a bound
        # for the whole asymptotic interval.
        self.assertAlmostEqual(float(bessel_j0(57.5)), decimal_j0(57.5),
                               delta=8e-15)


class AcousticVentReviewTests(unittest.TestCase):
    def test_dc_cutoff_and_normal_frequency_response(self):
        np.testing.assert_allclose(acoustic_vent_highpass([0, 25, 50], 25),
                                   [0, 0.5 + 0.5j, 0.8 + 0.4j],
                                   rtol=0, atol=2e-16)

    def test_finite_extreme_ratios_have_stable_limits(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            high = acoustic_vent_highpass([1e308], 1e-308)[0]
            low = acoustic_vent_highpass([1e-308], 1e308)[0]
        self.assertEqual(high, 1)
        self.assertEqual(low, 0)


class SeparationMetricReviewTests(unittest.TestCase):
    def test_ci_sdr_scores_the_full_convolution_tail(self):
        reference = np.array([1., 0., -1., 0.])
        estimate = np.array([1., -1., 1., -1.])
        result = ci_sdr_toeplitz(estimate, reference, 4)
        expected_filter = np.array([1., -1., 2., -2.]) / 3
        np.testing.assert_allclose(result["filter_taps"], expected_filter, atol=1e-15)
        np.testing.assert_allclose(result["target_projection"],
                                   np.convolve(reference, expected_filter), atol=1e-15)
        self.assertAlmostEqual(result["target_energy"], 4/3)
        self.assertAlmostEqual(result["residual_energy"], 8/3)
        self.assertAlmostEqual(result["ci_sdr_db"], -10 * math.log10(2))
        self.assertEqual(result["evaluation_samples"], 7)
        # The old cropped fit remains an explicitly different objective.
        cropped = causal_fir_toeplitz(reference, 4, full_convolution=False)
        np.testing.assert_array_equal(cropped @ [1, -1, 2, -2], estimate)

    def test_full_convolution_matrix_matches_independent_direct_convolution(self):
        reference, taps = np.array([2., -1., 3.]), np.array([1., 4., -2.])
        actual = causal_fir_toeplitz(reference, 3)
        self.assertEqual(actual.shape, (5, 3))
        np.testing.assert_array_equal(actual @ taps, np.convolve(reference, taps))

    def test_gate_keeps_large_finite_ratio_and_reports_energy_overflow(self):
        slots = np.array([[1e308, -1e308], [1e307, -1e307]])
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            result = single_speaker_bypass_gate(slots)
        self.assertAlmostEqual(result["secondary_to_primary_ratio"], 0.01)
        self.assertAlmostEqual(result["secondary_to_primary_db"], -20)
        self.assertTrue(result["bypass"])
        self.assertEqual(result["energy_representation"], ["overflow", "overflow"])
        self.assertIsNone(result["primary_energy"])
        np.testing.assert_allclose(result["output_waveforms"][0], [1.1e308, -1.1e308])

    def test_gate_preserves_log_ratio_when_linear_ratio_underflows(self):
        result = single_speaker_bypass_gate([[1e308, 0.], [1e-308, 0.]])
        self.assertTrue(result["ratio_underflow"])
        self.assertAlmostEqual(result["secondary_to_primary_db"], -12320)
        self.assertEqual(result["energy_representation"], ["overflow", "underflow"])
        self.assertTrue(result["bypass"])

    def test_gate_rejects_an_unrepresentable_merged_waveform(self):
        top = np.finfo(float).max
        with self.assertRaisesRegex(ValueError, "merged slot"):
            single_speaker_bypass_gate([[top], [top * 0.01]])


class TrackingReviewTests(unittest.TestCase):
    def test_psd_acceptance_is_scale_consistent(self):
        from codes.chapters.ch09.core.echo_av_fusion import _checked_psd_2x2
        for scale in (1e-200, 1e-16, 1., 1e200):
            with self.assertRaises(ValueError):
                _checked_psd_2x2(np.diag([-scale, scale]))
            np.testing.assert_array_equal(_checked_psd_2x2(np.diag([0., scale])),
                                          np.diag([0., scale]))
        with self.assertRaises(ValueError):
            _checked_psd_2x2(np.array([[1, 1j], [-1j, 1]]))
        maximum = np.finfo(float).max
        np.testing.assert_array_equal(_checked_psd_2x2(np.eye(2)*maximum),
                                      np.eye(2)*maximum)

    def test_ospa_tiny_positive_distance_is_preserved_in_both_entries(self):
        from codes.chapters.ch09.core.echo_av_fusion import angular_ospa_distance
        from codes.chapters.ch09.chapter09_experiments import small_set_distances
        for delta in (1e-200, 1., 9.):
            expected = delta
            self.assertEqual(angular_ospa_distance([0], [delta], order=2), expected)
            self.assertEqual(small_set_distances([0], [delta], order=2)['ospa'], expected)
        self.assertEqual(angular_ospa_distance([], [], order=2), 0)
        self.assertEqual(angular_ospa_distance([], [30], order=2), 10)
        self.assertEqual(angular_ospa_distance([179], [-179], order=2), 2)

    def test_lifecycle_rejects_skipped_frames_and_declares_hop(self):
        from codes.chapters.ch09.core.lifecycle import teaching_lifecycle
        def frames(starts):
            return {'start_sample': starts, 'observation_valid': [True]*len(starts),
                    'state_time_s': [(2*x+511)/32000 for x in starts],
                    'available_time_s': [2*(x+512)/32000 for x in starts]}
        with self.assertRaisesRegex(ValueError, 'hop_samples'):
            teaching_lifecycle(frames([0,1600,3200]))
        result = teaching_lifecycle(frames([0,160,320]))
        self.assertEqual(result['parameters']['hop_samples'], 160)
        self.assertEqual(result['events'][-1]['event'], 'confirmed')
        result = teaching_lifecycle(frames([0,80,160]), hop_samples=80)
        self.assertEqual(result['events'][-1]['frame'], 2)

    def test_mht_symmetric_branch_predictions_remain_tied(self):
        branches = (([41.,39.], [10.,-10.]), ([39.,41.], [10.,-10.]))
        squared_errors = [float(np.sum((np.array(pos)+vel-[50.,30.])**2))
                          for pos, vel in branches]
        self.assertEqual(squared_errors, [2., 2.])

    def test_miss_penalty_and_maximum_cardinality_are_distinct_objectives(self):
        # Two gated matches cost 16; one zero-cost match and one miss cost 9.
        full_matching = 8+8
        optional_missed_detection = 0+9
        self.assertLess(optional_missed_detection, full_matching)


class AecWpePremiseReviewTests(unittest.TestCase):
    def test_projection_trace_does_not_set_every_vector_energy(self):
        projection = np.diag([1., 0.])
        self.assertEqual(np.trace(projection)/2, .5)
        self.assertEqual(np.sum((projection @ [1.,0.])**2), 1.)
        self.assertEqual(np.sum((projection @ [0.,1.])**2), 0.)

    def test_floor_derivative_sign_not_global_increase(self):
        a = 2.
        first = lambda value: (value-a)/value**2
        self.assertGreater(first(3.), 0)
        self.assertGreater(first(5.), first(10.))
        candidates = np.array([3., 4., 6.])
        values = a/candidates + np.log(candidates)
        self.assertEqual(int(np.argmin(values)), 0)

    def test_fixed_loading_has_extra_identity_term(self):
        old = np.array([[2.,.5],[.5,1.]])
        q = np.array([1.,2.]); alpha, loading, power = .7, .4, 2.
        actual = alpha*old + np.outer(q,q)/power + loading*np.eye(2)
        adjusted = alpha*(old+loading*np.eye(2)) + np.outer(q,q)/power
        np.testing.assert_allclose(actual-adjusted, (1-alpha)*loading*np.eye(2), atol=1e-15)

    def test_common_gain_preserves_spatial_ratio_but_not_fixed_temporal_model(self):
        frame = np.array([1.+2j, 3.-1j])
        np.testing.assert_allclose((.3*frame)[1]/(.3*frame)[0], frame[1]/frame[0])
        self.assertNotEqual(.3*2., .1*2.)

    def test_shared_leading_delay_changes_mint_target_feasibility(self):
        matrix = np.array([[0.,0.],[1.,2.],[1.,-1.]])
        early = np.array([1.,0.,0.]); delayed = np.array([0.,1.,0.])
        early_fit = matrix @ np.linalg.lstsq(matrix, early, rcond=None)[0]
        delayed_fit = matrix @ np.linalg.lstsq(matrix, delayed, rcond=None)[0]
        self.assertEqual(float(np.linalg.norm(early_fit-early)), 1.)
        np.testing.assert_allclose(delayed_fit, delayed, atol=1e-15)

    def test_positive_power_need_not_have_positive_log_and_small_ridge_is_large_inverse(self):
        self.assertLess(math.log(.25), 0)
        self.assertEqual(1/1e-8, 1e8)

    def test_icc_report_names_listener_arrival_difference(self):
        from codes.chapters.ch06.core.nonlinear_and_stereo import stereo_decorrelation_and_loop_deadlines
        result = stereo_decorrelation_and_loop_deadlines()
        cases = result['loop_deadlines']['icc_itu_t_p1150']['cases']
        for case, expected in zip(cases, (10.,25.)):
            self.assertEqual(case['reinforced_first_arrival_ms']-case['direct_first_arrival_ms'], expected)
            self.assertEqual(case['arrival_difference_ms'], expected)
            self.assertEqual(case['meets_arrival_difference_limit'], expected<15)
            self.assertNotIn('itu_t_p1150_compliant', case)


class IndependentR3FollowupTests(unittest.TestCase):
    def test_raw_gram_underflow_is_explicit_without_losing_valid_projection(self):
        for scale, representation in ((1., 'finite'), (1e-150, 'finite'),
                                       (1e-200, 'positive_diagonal_underflow'),
                                       (1e200, 'overflow')):
            result = ci_sdr_toeplitz([1., -1.], scale*np.array([1., -1.]), 1)
            self.assertEqual(result['gram_matrix_representation'], representation)
            self.assertEqual(result['rank'], 1)
            self.assertEqual(result['ci_sdr_db'], 120.)
            np.testing.assert_allclose(result['target_projection'], [1.,-1.], atol=1e-15)
            if representation == 'finite':
                self.assertAlmostEqual(result['gram_matrix'][0,0]/(scale*scale), 2.)
            else:
                self.assertIsNone(result['gram_matrix'])

    def test_unrepresentable_positive_ci_sdr_energy_is_not_reported_as_zero(self):
        with self.assertRaisesRegex(ValueError, 'positive unscaled CI-SDR energy'):
            ci_sdr_toeplitz([1e-200,-1e-200], [1.,-1.], 1)

    def test_clear_remaining_assignment_uses_real_distance_before_exact_ties(self):
        from codes.chapters.ch09.core.echo_av_fusion import clear_mot_angular_summary
        first = 100.; second = np.nextafter(first, math.inf)
        result = clear_mot_angular_summary([{1:first,2:second}], [{1:second,2:first}])
        self.assertEqual(result['frames'][0]['matches'], {'1':2,'2':1})
        self.assertEqual(result['motp_deg'], 0.)
        tied = clear_mot_angular_summary([{1:0.,2:2.}], [{7:1.,8:1.}])
        self.assertEqual(tied['frames'][0]['matches'], {'1':7,'2':8})
        self.assertEqual(tied['motp_deg'], 1.)

    def test_principal_interval_preserves_adjacent_float_angles_and_set_members(self):
        from codes.chapters.ch09.core.tracking import wrap_angle
        from codes.chapters.ch09.core.echo_av_fusion import angular_ospa_distance, small_set_distances
        first = 100.; second = np.nextafter(first, math.inf)
        self.assertEqual(wrap_angle(first), first)
        self.assertEqual(wrap_angle(second), second)
        for order in (1,2):
            self.assertEqual(angular_ospa_distance([first], [second], order=order), second-first)
            self.assertEqual(small_set_distances([first], [second], order=order)['ospa'], second-first)
            self.assertEqual(small_set_distances([first,second], [first,second], order=order)['ospa'], 0.)
        np.testing.assert_array_equal(wrap_angle([-540,-180,180,540]), [-180]*4)
        for value in (np.finfo(float).max, -np.finfo(float).max):
            self.assertTrue(-180 <= wrap_angle(value) < 180)


class SpatialNoiseAndPhaseReviewTests(unittest.TestCase):
    def test_uca_radial_equalizer_requires_colored_noise_whitening(self):
        from codes.chapters.ch04.core.uca_cdr import uca_phase_mode_esprit, uca_steering_vector
        matrix = np.column_stack([uca_steering_vector(32, 1., angle) for angle in (-20.,45.)])
        signal = matrix @ np.diag([1.,.7]) @ matrix.conj().T
        for noise in (0.,1.,10.):
            estimates, _ = uca_phase_mode_esprit(signal+noise*np.eye(32), 1., n_sources=2)
            np.testing.assert_allclose(estimates, [-20.,45.], atol=1e-9)

    def test_dma_finite_extreme_input_returns_finite_weights_or_explicit_range_error(self):
        from codes.chapters.ch05.core.dma_kronecker import first_order_dma_weights
        for frequency, spacing, speed in ((1e308,.01,343.), (1e308,1e308,1.),
                                           (1e-308,1e-308,343.)):
            with warnings.catch_warnings():
                warnings.simplefilter('error', RuntimeWarning)
                try:
                    weights = first_order_dma_weights(frequency, spacing, 180., c=speed)
                except ValueError:
                    continue
            self.assertTrue(np.all(np.isfinite(weights)))


if __name__ == "__main__":
    unittest.main()
