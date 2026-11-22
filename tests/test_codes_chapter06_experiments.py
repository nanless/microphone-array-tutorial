"""Hand fractions and independent physical/metric anchors for E06-12..20,26..27,29..32."""
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
            {f'E06-{n:02d}' for n in range(12, 21)} | {'E06-26', 'E06-27', 'E06-29', 'E06-30', 'E06-31', 'E06-32'},
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

    def test_e06_31_pbfdkf_double_talk_variance_and_forgetting(self):
        r = self.results['E06-31']
        json.dumps(r, allow_nan=False)
        reg = r['regimes']
        # (a) Clean far-end single-talk: K = 4/5, W1 = 24/25, err1 = 1/25, P1 = 1/125.
        st = reg['clean_far_end_single_talk']
        self.assertAlmostEqual(st['kalman_gain'], float(F(4, 5)), places=15)
        self.assertAlmostEqual(st['posterior_weight'], float(F(24, 25)), places=15)
        self.assertAlmostEqual(st['posterior_error'], float(F(1, 25)), places=15)
        self.assertAlmostEqual(st['posterior_variance'], float(F(1, 125)), places=15)
        self.assertAlmostEqual(st['error_change_db'], 20.0 * math.log10(0.04 / 0.2), places=12)
        self.assertAlmostEqual(st['normalized_misalignment_db'], 20.0 * math.log10(0.04), places=12)
        # (b) Unguarded double-talk: K = 4/5, W1 = 64/25, err1 = -39/25 (+17.84 dB blowup).
        ug = reg['unguarded_double_talk']
        self.assertAlmostEqual(ug['residual_e'], 2.2, places=15)
        self.assertAlmostEqual(ug['kalman_gain'], float(F(4, 5)), places=15)
        self.assertAlmostEqual(ug['posterior_weight'], float(F(64, 25)), places=15)
        self.assertAlmostEqual(ug['posterior_error'], float(F(-39, 25)), places=15)
        self.assertAlmostEqual(ug['error_change_db'], 20.0 * math.log10(1.56 / 0.2), places=12)
        # (c) Adapted PBFDKF (Psi_e = |e|^2 = 4.84): K = 1/122, dW = 11/610, W1 = 499/610, err1 = 111/610.
        ad = reg['adapted_pbfdkf_double_talk']
        self.assertAlmostEqual(ad['psi_e'], 4.84, places=14)
        self.assertAlmostEqual(ad['kalman_gain'], float(F(1, 122)), places=14)
        self.assertAlmostEqual(ad['weight_increment'], float(F(11, 610)), places=14)
        self.assertAlmostEqual(ad['posterior_weight'], float(F(499, 610)), places=14)
        self.assertAlmostEqual(ad['posterior_error'], float(F(111, 610)), places=14)
        self.assertLess(ad['abs_posterior_error'], 0.2)
        # ML-subtracted variance (Psi_e = 4.80) and opposite-polarity burst (s = -2.0).
        ml = reg['adapted_pbfdkf_ml_variance']
        self.assertAlmostEqual(ml['kalman_gain'], float(F(1, 121)), places=15)
        self.assertAlmostEqual(ml['posterior_error'], float(F(2, 11)), places=15)
        opp = reg['adapted_pbfdkf_opposite_polarity']
        self.assertAlmostEqual(opp['kalman_gain'], float(F(1, 82)), places=15)
        self.assertAlmostEqual(opp['weight_increment'], float(F(-9, 410)), places=15)
        # Steady-state Riccati fixed point: A = 0.999 vs A = 0.99.
        ss = r['steady_state_forgetting']
        c999 = ss['cases']['A_0_999']
        c99 = ss['cases']['A_0_99']
        for case in (c999, c99):
            a_val = case['state_transition_A']
            p_m = case['prior_variance_ss']
            p_p = case['posterior_variance_ss']
            k_ss = case['kalman_gain_ss']
            psi_d = case['process_variance_psi_delta']
            self.assertAlmostEqual((a_val ** 2) * p_p + psi_d, p_m, places=14)
            self.assertAlmostEqual(p_p, k_ss * ss['psi_e'], places=14)
        self.assertGreater(c99['kalman_gain_ss'], 2.0 * c999['kalman_gain_ss'])
        self.assertGreater(ss['excess_floor_increase_db'], 3.0)

    def test_e06_32_res_wiener_gain_floor_and_cng_matching(self):
        r = self.results['E06-32']
        json.dumps(r, allow_nan=False)
        self.assertAlmostEqual(r['g_min'], 0.1, places=15)
        self.assertAlmostEqual(r['target_background_noise_psd'], 0.01, places=15)
        # (a) Far-end single-talk: unclamped 0, clamped 0.1 (-20 dB), noise 1e-4 (-40 dB), CNG 0.0099.
        fe = r['frames']['far_end_single_talk']
        self.assertEqual(fe['g_unclamped'], 0.0)
        self.assertAlmostEqual(fe['g_res'], 0.1, places=15)
        self.assertTrue(fe['clamped_by_floor'])
        self.assertAlmostEqual(fe['echo_suppression_db'], 20.0, places=12)
        self.assertAlmostEqual(fe['output_echo_psd'], 0.0025, places=15)
        self.assertAlmostEqual(fe['output_noise_psd_without_cng'], 1e-4, places=15)
        self.assertAlmostEqual(fe['noise_pumping_db_without_cng'], -20.0, places=12)
        self.assertAlmostEqual(fe['cng_psd'], 0.0099, places=15)
        self.assertAlmostEqual(fe['output_noise_psd_with_cng'], 0.01, places=15)
        self.assertAlmostEqual(fe['noise_pumping_db_with_cng'], 0.0, places=12)
        # (b) Double-talk: G = 50/63, noise_no_cng = 25/3969, cng = 1469/396900 -> exact 0.01.
        dt = r['frames']['double_talk']
        self.assertFalse(dt['clamped_by_floor'])
        self.assertAlmostEqual(dt['g_res'], float(F(50, 63)), places=15)
        self.assertAlmostEqual(dt['speech_attenuation_db'], 20.0 * math.log10(50.0 / 63.0), places=12)
        self.assertAlmostEqual(dt['echo_suppression_db'], -20.0 * math.log10(50.0 / 63.0), places=12)
        self.assertAlmostEqual(dt['output_noise_psd_without_cng'], float(F(25, 3969)), places=15)
        self.assertAlmostEqual(dt['cng_psd'], float(F(1469, 396900)), places=15)
        self.assertAlmostEqual(dt['output_noise_psd_with_cng'], 0.01, places=15)
        self.assertAlmostEqual(dt['noise_pumping_db_with_cng'], 0.0, places=12)


if __name__ == '__main__':
    unittest.main()
