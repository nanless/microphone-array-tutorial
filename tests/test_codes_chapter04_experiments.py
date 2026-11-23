"""Independent closed-form anchors and numerical failure checks for chapter 4."""
import cmath
import json
import math
import unittest
from scripts import run_portable_ci as portable_ci
import warnings
import numpy as np
from codes.chapters.ch04 import chapter04_experiments as ex
from codes.chapters.ch04.core.doa import esprit_ula, gcc_phat
from codes.chapters.ch04.core.uca_cdr import (
    FIRST_BESSEL_J0_ROOT,
    FIRST_BESSEL_J1_ROOT,
    bessel_jn,
    cdr_nodir,
    uca_phase_mode_esprit,
    uca_steering_vector,
)
from codes.chapters.ch00.core.audio_samples import doa_ambiguity_case, prepare_exports, read_pcm16


class Chapter04ExperimentsTests(unittest.TestCase):
    def test_catalog_keys_and_strict_json(self):
        if not portable_ci.is_original_generation_host():
            self.skipTest("complete published exercise replay requires the original generation host")
        results = ex.run_exercises()
        self.assertEqual(set(results), {f'E04-{i:02d}' for i in range(12, 25)})
        json.dumps(results, allow_nan=False)

    def test_nyquist_recovers_original_grid_and_half_sample(self):
        r = ex.nyquist_interpolation()
        y = np.array(r['correlation'])
        np.testing.assert_allclose(y[::2], [0, 0, 0, 0, 1, 0, 0, 0], atol=2e-15)
        # Dirichlet endpoint splitting: cot(pi/16)/8 at half a sample.
        self.assertAlmostEqual(y[9], 1/(8*np.tan(np.pi/16)), places=14)
        self.assertLess(r['max_error'], 2e-15)

    def test_srp_two_frequency_scalar_cosines(self):
        r = ex.srp_frequency_and_grid()
        np.testing.assert_allclose(r['far_field_scores'], [-.5, np.sqrt(2)/4, 1], atol=1e-14)
        np.testing.assert_allclose(r['grid_peak_m'], [5.567415730337078, 3.838805970149254])
        self.assertAlmostEqual(r['grid_error_m'], .0777867856136663)
        self.assertGreater(r['grid_peak_score'], r['nearest_score'])
        p = r['grid_peak_m']; source = [5.5, 3.8]
        mics = [(0, 0), (4, 0), (0, 3), (4, 3)]
        delta = [math.dist(p, m)-math.dist(source, m) for m in mics]
        reference = sum(math.exp(-.5*((delta[i]-delta[j])*6000/343)**2)
                        for i in range(4) for j in range(i+1, 4))/6
        self.assertAlmostEqual(r['grid_peak_score'], reference, places=14)

    def test_phat_order_and_epsilon_boundaries(self):
        r = ex.phat_order_and_floors()
        self.assertEqual((r['framewise_then_average'], r['average_then_normalize']), (0, -1))
        np.testing.assert_allclose(r['add'], [100/101, .5, 1/11, 0])
        np.testing.assert_allclose(r['floor'], [1, 1, .1, 0])
        np.testing.assert_array_equal(r['gate'], [1, 0, 0, 0])

    def test_esprit_common_basis_and_principal_alias(self):
        r = ex.esprit_basis_and_alias()
        psi = np.array(r['mixed_basis_shift_real'])+1j*np.array(r['mixed_basis_shift_imag'])
        np.testing.assert_allclose(psi, [[-1j, -2j], [0, 1j]], atol=1e-14)
        np.testing.assert_allclose(r['estimated_angles_deg'], [-30, 30], atol=1e-12)
        self.assertAlmostEqual(r['principal_alias_deg'], np.rad2deg(np.arcsin(np.sin(np.deg2rad(40))-1)))

    def test_esprit_rejects_uninformative_nonpsd_rank_and_zero_mode(self):
        zero_mode = np.outer([1., 0., 1.], [1., 0., 1.])
        for matrix in [np.zeros((3, 3)), np.diag([-2., -1., 3.]),
                       np.diag([0., 0., 1.]), np.eye(3), zero_mode]:
            with self.subTest(matrix=matrix), self.assertRaises(np.linalg.LinAlgError):
                esprit_ula(matrix, source_count=1, spacing_m=.03, frequency_hz=1000)
        with self.assertRaises(ValueError):
            esprit_ula(np.eye(3), source_count=True, spacing_m=.03, frequency_hz=1000)

    def test_esprit_rejects_overselected_signal_dimension(self):
        a = np.exp(1j*.2*np.arange(4))
        for floor in [0., .1]:
            with self.subTest(floor=floor), self.assertRaises(np.linalg.LinAlgError):
                esprit_ula(np.outer(a, a.conj())+floor*np.eye(4), source_count=2,
                           spacing_m=.03, frequency_hz=1000)

    def test_esprit_extreme_positive_covariance_scaling(self):
        a = np.array([1., 1j, -1., -1j])
        covariance = np.outer(a, a.conj())+.1*np.eye(4)
        for factor in [1e-250, 1., 1e250]:
            with self.subTest(factor=factor):
                actual = esprit_ula(factor*covariance, source_count=1, spacing_m=.1715, frequency_hz=1000)
                np.testing.assert_allclose(actual, [np.pi/6], atol=1e-14)

    def test_gcc_positive_gain_is_not_polarity_invariance(self):
        reference = np.zeros(16); reference[5] = 1
        delayed = np.zeros(16); delayed[7] = 1
        original = gcc_phat(delayed, reference, 16000)
        scaled = gcc_phat(.001*delayed, 1000*reference, 16000)
        self.assertAlmostEqual(original[0], 2/16000)
        self.assertAlmostEqual(scaled[0], original[0])
        self.assertNotAlmostEqual(gcc_phat(-delayed, reference, 16000)[0], original[0])

    def test_gauss_newton_closed_form_normal_equations(self):
        combined = ex.gauss_newton_and_gls()
        r = combined['gauss_newton_step']
        z = np.array([np.sqrt(.85)-np.sqrt(.65), .5-np.sqrt(.65), np.sqrt(.45)-np.sqrt(.65)])
        j = -np.sqrt(2)*np.array([[1., 0.], [0., 1.], [1., 1.]])
        expected = np.array([[2., -1.], [-1., 2.]])/6 @ j.T @ z
        np.testing.assert_allclose(r['observation_range_difference_m'], z)
        np.testing.assert_allclose(r['jacobian'], j)
        np.testing.assert_allclose(r['step_m'], expected)
        self.assertAlmostEqual(r['residual_after_m'], .013927510559068, places=13)
        gls = combined['correlated_tdoa_gls']['methods']
        np.testing.assert_allclose(gls['gls']['step_m'], [-.03/math.sqrt(2), 0.], atol=1e-15)

    def test_crlb_exact_sine_square_sum(self):
        r = ex.local_delay_crlb()
        self.assertAlmostEqual(r['fisher_information_per_s2']/(400_000_000*np.pi**2), 4)
        self.assertAlmostEqual(r['std_bound_us'], 25/np.pi, places=12)
        self.assertAlmostEqual(r['variance_bound_s2']*1e12, 625/np.pi**2, places=12)

    def test_audio_actual_pcm_steady_delay_and_quantization(self):
        case = doa_ambiguity_case()
        files, groups = prepare_exports({'doa_ambiguity': case})
        self.assertEqual(groups['doa_ambiguity']['common_export_gain'], 1)
        for filename, (blob, metadata) in files.items():
            fs, pcm = read_pcm16(blob)
            self.assertEqual((fs, pcm.shape), (16000, (2, 32000)))
            np.testing.assert_array_equal(pcm[0, 2:], pcm[1, :-2])
            np.testing.assert_array_equal(pcm[0, :2], [0, 0])
            self.assertLessEqual(metadata['quantization_max_abs_error'], .5/32768)

    def test_audio_scores_geometric_series_and_physical_aliases(self):
        r = ex.audio_phase_ambiguity()
        np.testing.assert_allclose(r['rows'][0]['pcm_scores'], [1, 1], atol=1e-14)
        a = -np.pi/64
        total = np.exp(1j*20*a)*(1-np.exp(1j*365*a))/(1-np.exp(1j*a))
        np.testing.assert_allclose(r['rows'][1]['pcm_scores'], [total.real/365, 1], atol=1e-14)
        self.assertTrue(all(abs(x) < r['physical_max_delay_samples'] for x in [-6, 2]))
        for row in r['rows']:
            np.testing.assert_allclose(row['pcm_scores'], row['float_scores'], atol=1e-14)

    def test_bessel_jn_series_parity_recurrence_and_type_guards(self):
        # Independent 9-term alternating series at x = 1:
        j0_hand = sum(((-1)**k) * (0.25**k) / (math.factorial(k)**2) for k in range(9))
        j1_hand = sum(((-1)**k) * (0.5**(1 + 2*k)) / (math.factorial(k)*math.factorial(k + 1)) for k in range(9))
        self.assertAlmostEqual(float(bessel_jn(0, 0.0)), 1.0, places=15)
        self.assertAlmostEqual(float(bessel_jn(1, 0.0)), 0.0, places=15)
        self.assertAlmostEqual(float(bessel_jn(0, 1.0)), j0_hand, places=15)
        self.assertAlmostEqual(float(bessel_jn(1, 1.0)), j1_hand, places=15)
        self.assertAlmostEqual(float(bessel_jn(-1, 1.0)), -j1_hand, places=15)
        self.assertAlmostEqual(float(bessel_jn(-2, 1.0)), float(bessel_jn(2, 1.0)), places=15)
        # Recurrence J_{n-1}(x) + J_{n+1}(x) == (2n/x) * J_n(x) across series (1.5) and quadrature (15.0):
        for x_val in (1.5, 15.0):
            lhs = float(bessel_jn(1, x_val) + bessel_jn(3, x_val))
            rhs = (4.0 / x_val) * float(bessel_jn(2, x_val))
            self.assertAlmostEqual(lhs, rhs, places=13)
        self.assertLess(abs(float(bessel_jn(0, FIRST_BESSEL_J0_ROOT))), 1e-14)
        self.assertLess(abs(float(bessel_jn(1, FIRST_BESSEL_J1_ROOT))), 1e-14)
        # High-order small-magnitude series convergence (regression for premature 1e-17*max(1, acc) break):
        j45_ref = 3.946558989987071e-22
        j45_val = float(bessel_jn(45, 12.0))
        self.assertLess(abs(j45_val - j45_ref) / j45_ref, 1e-12)
        j45_arr = bessel_jn(45, np.array([0.0, 12.0]))
        self.assertEqual(float(j45_arr[0]), 0.0)
        self.assertLess(abs(float(j45_arr[1]) - j45_ref) / j45_ref, 1e-12)
        for bad_n in (True, 1.5, '0'):
            with self.assertRaises(ValueError):
                bessel_jn(bad_n, 1.0)
        for bad_x in (True, 1.0 + 0.0j, float('nan'), float('inf')):
            with self.assertRaises(ValueError):
                bessel_jn(0, bad_x)

    def test_e04_22_uca_phase_mode_esprit_aliasing_and_bessel_zeros(self):
        r = ex.exercise_04_22()
        self.assertAlmostEqual(r['j0_1'], 0.7651976865579666, places=15)
        self.assertAlmostEqual(r['j1_1'], 0.44005058574493355, places=15)
        self.assertLess(r['abs_bias_m8_deg'], 3e-4)
        # Independent 4-point DFT derivation for M = 4 circular harmonic aliasing at phi = pi/3:
        phi = math.pi / 3.0
        psi4 = [0.5 * math.pi * q for q in range(4)]
        a4 = [cmath.exp(1j * math.cos(phi - p)) for p in psi4]
        w_rad = [-1j / r['j1_1'], 1.0 / r['j0_1'], -1j / r['j1_1']]
        tilde4 = [
            w_rad[idx] * sum(cmath.exp(-1j * m * psi4[q]) * a4[q] for q in range(4)) / 4.0
            for idx, m in enumerate((-1, 0, 1))
        ]
        mu4 = (tilde4[0].conjugate()*tilde4[1] + tilde4[1].conjugate()*tilde4[2]) / (abs(tilde4[0])**2 + abs(tilde4[1])**2)
        expected_m4_deg = math.degrees(0.5 * math.pi + cmath.phase(mu4))
        self.assertAlmostEqual(r['theta_hat_m4_deg'], expected_m4_deg, places=12)
        self.assertGreater(r['abs_bias_m4_deg'], 2.0)
        self.assertAlmostEqual(r['f_zero_m0_hz'], 343.0 * FIRST_BESSEL_J0_ROOT / (0.08 * math.pi), places=10)
        self.assertAlmostEqual(r['f_zero_m1_hz'], 343.0 * FIRST_BESSEL_J1_ROOT / (0.08 * math.pi), places=10)
        self.assertTrue(r['bessel_singularity_rejected'])
        # Full 360-deg azimuth quadrant coverage (-135 deg, +140 deg) on M = 8:
        for target_az in (-135.0, 140.0):
            vec = uca_steering_vector(8, 1.0, target_az)
            est, _ = uca_phase_mode_esprit(np.outer(vec, vec.conj()) + 1e-4 * np.eye(8), 1.0, n_sources=1)
            self.assertAlmostEqual(float(est[0]), target_az, delta=1e-3)
        # Subnormal covariance scaling (1e-315) must avoid complex division overflow/NaN:
        with warnings.catch_warnings():
            warnings.simplefilter('error', RuntimeWarning)
            vec30 = uca_steering_vector(8, 1.0, 30.0)
            cov_sub = 1e-315 * (np.outer(vec30, vec30.conj()) + 1e-3 * np.eye(8))
            est_sub, _ = uca_phase_mode_esprit(cov_sub, 1.0, n_sources=1)
            self.assertAlmostEqual(float(est_sub[0]), 30.0, delta=1e-3)
        with self.assertRaises(np.linalg.LinAlgError):
            uca_phase_mode_esprit(np.eye(8), 1.0, n_sources=1)

    def test_e04_23_cdr_nodir_exact_rational_and_regimes(self):
        r = ex.exercise_04_23()
        # Hand-derived exact rationals: gamma_x = 1/6 + 2j/3, |gamma_x|^2 = 17/36, CDR = 2, w = 2/3.
        self.assertAlmostEqual(r['observed_gamma_x_real'], 1.0 / 6.0, places=15)
        self.assertAlmostEqual(r['observed_gamma_x_imag'], 2.0 / 3.0, places=15)
        self.assertAlmostEqual(r['observed_abs_gamma_x_sq'], 17.0 / 36.0, places=15)
        self.assertAlmostEqual(r['cdr_est'], 2.0, places=14)
        self.assertAlmostEqual(r['tf_weight'], 2.0 / 3.0, places=14)
        self.assertEqual(r['cdr_pure_diffuse'], 0.0)
        self.assertTrue(r['cdr_pure_anechoic_is_infinite'])
        self.assertEqual(r['tf_weight_pure_anechoic'], 1.0)
        self.assertTrue(r['dc_degeneracy_rejected'])
        with self.assertRaises(ValueError):
            cdr_nodir(1.01 + 0.0j, 0.5)
        with self.assertRaises(ValueError):
            cdr_nodir(0.5, 1.0)

    def test_e04_24_off_grid_l1_mismatch_and_seld_f20(self):
        r = ex.exercise_04_24()
        json.dumps(r, allow_nan=False)
        l1 = r['off_grid_l1_mismatch']
        self.assertEqual((l1['num_mics'], l1['grid_size'], l1['grid_step_deg']), (8, 25, 5.0))
        self.assertAlmostEqual(l1['mu'], 0.02, places=15)
        on = l1['on_grid']
        off = l1['off_grid']
        # On-grid exact analytical LASSO optimum: single atom at 10 deg with weight 1 - mu = 0.98.
        self.assertEqual(on['support_count'], 1)
        self.assertEqual(on['peak_angle_deg'], 10.0)
        self.assertAlmostEqual(on['peak_magnitude'], 0.98, places=12)
        self.assertAlmostEqual(on['peak_energy_fraction'], 1.0, places=14)
        self.assertAlmostEqual(on['weight_at_15_deg'], 0.0, places=14)
        self.assertEqual(on['sidelobe_atom_count'], 0)
        self.assertAlmostEqual(on['relative_residual_l2'], 0.02, places=12)
        self.assertLess(on['kkt_active_error'], 1e-12)
        self.assertLess(on['kkt_inactive_error'], 1e-12)
        # Off-grid 12.5 deg splits across 10 & 15 deg plus two spurious sidelobes (-5, +30 deg).
        self.assertEqual(off['support_count'], 4)
        self.assertEqual([a['angle_deg'] for a in off['active_atoms']], [-5.0, 10.0, 15.0, 30.0])
        self.assertEqual(off['peak_angle_deg'], 15.0)
        self.assertAlmostEqual(off['weight_at_10_deg'], 0.506397, places=5)
        self.assertAlmostEqual(off['weight_at_15_deg'], 0.515831, places=5)
        self.assertGreater(off['weight_at_15_deg'], off['weight_at_10_deg'])
        self.assertEqual(off['sidelobe_atom_count'], 2)
        self.assertAlmostEqual(off['adjacent_10_15_energy_fraction'], 0.999568, places=5)
        self.assertAlmostEqual(off['sidelobe_energy_fraction'], 0.000432, places=5)
        self.assertAlmostEqual(off['relative_residual_l2'], 0.034898, places=5)
        self.assertGreater(off['relative_residual_l2'], on['relative_residual_l2'])
        self.assertGreater(off['l1_norm'], 1.0)
        self.assertLess(off['kkt_active_error'], 1e-12)
        self.assertLess(off['kkt_inactive_error'], 1e-12)
        # Part 2: SELD location-dependent F_{20 deg} coupling.
        seld = r['seld_location_dependent_f1']
        self.assertEqual(seld['doa_threshold_deg'], 20.0)
        before = seld['before_refinement']
        after = seld['after_refinement']
        self.assertEqual((before['tp_20'], before['fp_20'], before['fn_20']), (3, 1, 1))
        self.assertAlmostEqual(before['f_20_percent'], 75.0, places=12)
        self.assertAlmostEqual(before['lr_cd_percent'], 100.0, places=12)
        self.assertAlmostEqual(before['le_cd_deg'], 15.0, places=12)
        self.assertEqual((after['tp_20'], after['fp_20'], after['fn_20']), (4, 0, 0))
        self.assertAlmostEqual(after['f_20_percent'], 100.0, places=12)
        self.assertAlmostEqual(after['lr_cd_percent'], 100.0, places=12)
        self.assertAlmostEqual(after['le_cd_deg'], 13.5, places=12)
        self.assertAlmostEqual(seld['f_20_gain_percentage_points'], 25.0, places=12)
        self.assertAlmostEqual(seld['le_cd_reduction_deg'], 1.5, places=12)


if __name__ == '__main__':
    unittest.main()
