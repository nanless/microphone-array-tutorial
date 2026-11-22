"""Hand fractions and independent physical/metric anchors for E06-12..20,26..27,29..30."""
from fractions import Fraction as F
import json
import math
import unittest
import numpy as np
from codes.chapters.ch06.chapter06_experiments import run_experiments
from codes.chapters.ch06.core import nonlinear_and_stereo as ns


class Chapter06Experiments(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results = run_experiments()

    def test_ids_and_finite_serialization(self):
        self.assertEqual(
            set(self.results),
            {f'E06-{n:02d}' for n in range(12, 21)} | {'E06-26', 'E06-27', 'E06-29', 'E06-30'},
        )
        json.dumps(self.results, allow_nan=False)

    def test_shared_normalization_independent_error_recurrence(self):
        r = self.results['E06-12']['normalization']
        self.assertEqual(r['joint']['prior_errors'], [1., -.5, .25])
        self.assertEqual(r['separate']['prior_errors'], [1., -2., 4.])
        self.assertEqual(r['joint']['next_error'], -.125)
        self.assertEqual(r['separate']['next_error'], -8.)

    def test_program_change_exposes_unidentified_difference(self):
        r = self.results['E06-13']
        self.assertEqual(r['correlation_eigenvalues'], [0., 2.])
        self.assertEqual(r['true_training_echo'], [3., -3.])
        self.assertEqual(r['alternative_training_echo'], [3., -3.])
        self.assertEqual(r['true_holdout_echo'], -1.)
        self.assertEqual(r['alternative_holdout_echo'], 0.)

    def test_delay_support_and_phase_from_quarter_period(self):
        c = self.results['E06-14']['cases']
        self.assertEqual([v['remaining_path_delay_samples'] for v in c], [2, 0, -1])
        self.assertEqual([v['representable_by_causal_fir'] for v in c], [True, True, False])
        r = self.results['E06-15']
        self.assertEqual(r['uncompensated_delay_samples'], [0, 4, 8])
        # 1 kHz has 16 samples per cycle; 4 and 8 samples are quarter/half cycle.
        np.testing.assert_allclose(r['residual_to_echo_power_ratio'], [0, 2, 4], atol=1e-15)

    def test_dropout_tail_uses_delayed_missing_samples(self):
        r = self.results['E06-16']
        self.assertEqual(r['echo'], [1, 2.5, 4, 5.5, 7, 8.5])
        self.assertEqual(r['prior_residual'], [0, 0, 3, 5.5, 2, 0])
        self.assertEqual(r['first_fully_recovered_output_sample'], 5)

    def test_prior_is_the_scored_output(self):
        r = self.results['E06-17']
        self.assertEqual((r['gain'], r['posterior_weight'], r['prior_error'], r['posterior_error']),
                         (.5, 1., 2., 1.))
        self.assertAlmostEqual(r['same_observation_fit_power_difference_db'], 10 * math.log10(4))

    def test_aggregate_score_is_not_mean_db(self):
        r = self.results['E06-18']
        self.assertEqual(r['arithmetic_mean_block_db'], 10.)
        self.assertAlmostEqual(r['combined_energy_ratio_db'], 10 * math.log10(float(F(10100, 10001))))

    def test_normalized_gradient_has_different_expectation(self):
        r = self.results['E06-19']
        self.assertEqual(r['cross_second_moment'], 0.)
        self.assertEqual(r['mean_near_over_reference'], float(F(3, 8)))
        self.assertAlmostEqual(r['batch_least_squares_path'], float(F(4, 5)))
        self.assertAlmostEqual(r['expected_nlms_equilibrium'], float(F(47, 40)))

    def test_increment_error_pythagorean_identity(self):
        r = self.results['E06-20']
        self.assertAlmostEqual(r['projection_gain'], float(F(4, 5)))
        self.assertAlmostEqual(r['target_remainder_inner_product'], 0.)
        self.assertAlmostEqual(r['gain_error_term'], float(F(1, 25)))
        self.assertAlmostEqual(r['orthogonal_error_term'], float(F(9, 100)))
        self.assertAlmostEqual(r['total_relative_squared_error'], float(F(13, 100)))

    def test_hammerstein_legendre_conditioning_and_decoupling(self):
        r = self.results['E06-29']
        expected_raw = [[float(F(1, p + q + 1)) for q in (1, 3, 5)] for p in (1, 3, 5)]
        np.testing.assert_allclose(r['raw_power_gram'], expected_raw, rtol=1e-15, atol=1e-15)
        self.assertAlmostEqual(r['raw_power_determinant'], float(F(256, 22920975)), places=18)
        self.assertAlmostEqual(
            r['raw_power_condition_number_2branch'],
            (25.0 + math.sqrt(541.0)) / (25.0 - math.sqrt(541.0)),
            places=12,
        )
        raw_eigs = r['raw_power_eigenvalues_ascending']
        for mu in raw_eigs:
            poly = (
                mu ** 3
                - float(F(131, 231)) * (mu ** 2)
                + float(F(19816, 1091475)) * mu
                - float(F(256, 22920975))
            )
            self.assertAlmostEqual(poly, 0.0, places=14)
        self.assertAlmostEqual(r['raw_power_condition_number_3branch'], 849.5758546, places=5)

        self.assertAlmostEqual(r['unit_power_determinant'], float(F(256, 99225)), places=15)
        self.assertAlmostEqual(
            r['unit_power_condition_number_2branch'],
            (5.0 + math.sqrt(21.0)) / (5.0 - math.sqrt(21.0)),
            places=12,
        )
        for lam in r['unit_power_eigenvalues_ascending']:
            poly = (
                lam ** 3
                - 3.0 * (lam ** 2)
                + float(F(53176, 99225)) * lam
                - float(F(256, 99225))
            )
            self.assertAlmostEqual(poly, 0.0, places=14)
        self.assertAlmostEqual(r['unit_power_condition_number_3branch'], 567.4541521, places=5)
        self.assertAlmostEqual(r['legendre_condition_number_3branch'], float(F(11, 3)), places=14)
        self.assertEqual(r['orthonormal_legendre_condition_number_3branch'], 1.0)

        np.testing.assert_allclose(
            r['legendre_coefficients'],
            [float(F(97, 100)), float(F(2, 25)), float(F(2, 25))],
            rtol=1e-15,
            atol=1e-15,
        )
        self.assertAlmostEqual(r['best_linear_truncation_gain'], float(F(97, 100)), places=15)
        self.assertTrue(r['linear_gain_matches_first_legendre_weight'])
        self.assertAlmostEqual(r['linear_explained_variance'], float(F(9409, 30000)), places=15)
        self.assertAlmostEqual(r['residual_nonlinear_variance'], float(F(72, 48125)), places=15)
        self.assertAlmostEqual(r['total_output_variance'], float(F(727949, 2310000)), places=15)
        self.assertLessEqual(r['max_pointwise_basis_difference'], 1e-15)

    def test_stereo_tv_allpass_deadlines_and_switching_contracts(self):
        r = self.results['E06-30']
        lti = r['stereo_decorrelation']['lti_allpass']
        np.testing.assert_allclose(
            lti['x2'],
            [float(F(1, 2)), float(F(3, 4)), float(F(-3, 8)), float(F(11, 16)), float(F(21, 32))],
            atol=1e-15,
        )
        self.assertAlmostEqual(lti['determinant'], 0.0, places=15)
        self.assertEqual(lti['rank'], 3)
        np.testing.assert_allclose(lti['null_product'], [0.0, 0.0, 0.0, 0.0], atol=1e-15)

        tv = r['stereo_decorrelation']['time_varying_allpass']
        np.testing.assert_allclose(tv['x2'], [0.0, 1.0, 0.0, 0.5, 0.75], atol=1e-15)
        self.assertAlmostEqual(tv['determinant'], 0.25, places=15)
        self.assertEqual(tv['rank'], 4)
        np.testing.assert_allclose(tv['echo_snapshots'], [0.8, 0.6, 1.0, 0.55], atol=1e-15)
        np.testing.assert_allclose(tv['recovered_stereo_path'], [0.8, 0.2, 0.4, 0.1], atol=1e-15)

        anc = r['loop_deadlines']['anc']
        self.assertAlmostEqual(anc['max_electronic_delay_ms'], float(F(100, 343)), places=12)
        self.assertAlmostEqual(anc['max_electronic_delay_samples'], float(F(1600, 343)), places=11)
        self.assertAlmostEqual(anc['block_to_budget_ratio'], 27.44, places=12)
        self.assertFalse(anc['broadband_causal_anc_feasible'])
        self.assertTrue(r['loop_deadlines']['digital_aec']['ram_reference_alignment_feasible'])

        icc = r['loop_deadlines']['icc_itu_t_p1150']
        self.assertAlmostEqual(icc['comb_peak_to_notch_swing_db'], 20.0 * math.log10(9.0), places=12)
        self.assertAlmostEqual(icc['howling_gain_margin_db'], 20.0 * math.log10(2.0), places=12)
        c10, c25 = icc['cases']
        self.assertEqual(
            (c10['itu_t_p1150_compliant'], c10['comb_notch_spacing_hz'], c10['first_notch_hz'], c10['notches_in_0_to_1000_hz']),
            (True, 100.0, 50.0, 10),
        )
        self.assertEqual(
            (c25['itu_t_p1150_compliant'], c25['comb_notch_spacing_hz'], c25['first_notch_hz'], c25['notches_in_0_to_1000_hz']),
            (False, 40.0, 20.0, 25),
        )

        wins = r['dual_filter_switching']['windows']
        self.assertEqual(wins['window_A_path_change']['action'], 'copy_shadow_to_foreground')
        self.assertEqual(wins['window_B_double_talk_diverged']['action'], 'reset_shadow_from_foreground')
        self.assertEqual(wins['window_C_double_talk_overfit_blocked']['action'], 'hold')

        with self.assertRaises(ValueError):
            ns.first_order_allpass_decorrelate([1.0, 0.0], [1.0, 0.0])
        with self.assertRaises(ValueError):
            ns.dual_filter_switch_decision(True, 0.1, 0.9)


if __name__ == '__main__':
    unittest.main()
