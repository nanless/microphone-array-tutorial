"""Independent hand-derived unit tests and strict input-validation negative controls for Chapter 0."""
from __future__ import annotations

import decimal
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys
import unittest

import numpy as np

from codes.chapters.ch00 import chapter00_exercises as ex00
from codes.chapters.ch00.core import chapter00_math as cm
from codes.chapters.ch00.io_contracts import strict_json_loads

ROOT = Path(__file__).resolve().parents[1]


class Chapter00MathKernelTests(unittest.TestCase):
    def test_e01_geometry_and_fresnel_hand_anchors(self):
        res = cm.spherical_and_fresnel_geometry_example()
        cb = res["coordinate_bridge"]
        expected_u = np.array([math.sqrt(2.0) / 4.0, math.sqrt(6.0) / 4.0, math.sqrt(2.0) / 2.0])
        np.testing.assert_allclose(cb["unit_direction"], expected_u, rtol=0.0, atol=1e-15)
        np.testing.assert_allclose(cb["position_cartesian_m"], 1.5 * expected_u, rtol=0.0, atol=1e-15)
        self.assertAlmostEqual(cb["alpha_ccw_from_px_deg"], 60.0, places=13)
        self.assertAlmostEqual(cb["colatitude_deg"], 45.0, places=13)
        self.assertAlmostEqual(cb["sofa_azimuth_deg"], -30.0, places=13)
        self.assertAlmostEqual(cb["sofa_elevation_deg"], 45.0, places=13)
        expected_u_sofa = np.array([math.sqrt(6.0) / 4.0, -math.sqrt(2.0) / 4.0, math.sqrt(2.0) / 2.0])
        np.testing.assert_allclose(cb["ambix_sofa_unit_direction"], expected_u_sofa, rtol=0.0, atol=1e-15)
        np.testing.assert_allclose(cb["ambix_sofa_position_cartesian_m"], 1.5 * expected_u_sofa, rtol=0.0, atol=1e-15)
        self.assertLess(cb["ambix_sofa_consistency_error"], 1e-15)
        self.assertAlmostEqual(res["spherical_law_of_cosines"]["central_angle_deg"], 90.0, places=13)

        fr = res["fresnel_approximation"]
        self.assertAlmostEqual(fr["wavelength_m"], 0.10, places=15)
        self.assertAlmostEqual(fr["fraunhofer_distance_m"], 0.80, places=15)
        self.assertTrue(fr["in_nearfield_fresnel_zone"])
        self.assertAlmostEqual(fr["delta_r_exact_m"], math.sqrt(0.21) - 0.50, places=15)
        self.assertAlmostEqual(fr["delta_r_plane_wave_m"], -0.05, places=15)
        self.assertAlmostEqual(fr["delta_r_fresnel_m"], -0.0425, places=15)

    def test_e02_so3_noncommutativity_and_norm_preservation(self):
        res = cm.so3_rotation_matrix_composition_example()
        self.assertAlmostEqual(res["determinant"], 1.0, places=15)
        self.assertLess(res["orthogonality_fro_error"], 1e-15)
        np.testing.assert_allclose(res["ex_under_rz_then_ry"], [0.0, 0.0, -1.0], rtol=0.0, atol=1e-15)
        np.testing.assert_allclose(res["ex_under_ry_then_rz"], [0.0, 1.0, 0.0], rtol=0.0, atol=1e-15)
        self.assertAlmostEqual(res["commutator_frobenius_norm"], math.sqrt(6.0), places=14)
        np.testing.assert_allclose(res["rotated_vector"], [-2.0, 2.0, -1.0], rtol=0.0, atol=1e-15)
        self.assertAlmostEqual(res["rotated_norm"], 3.0, places=15)
        self.assertLess(res["inverse_recovery_error"], 1e-15)

    def test_e03_phasor_delay_and_db_budgets(self):
        res = cm.phasor_delay_and_decibel_budget_example()
        self.assertAlmostEqual(res["quadrature_phasor"].real, 1.0, places=15)
        self.assertAlmostEqual(res["quadrature_phasor"].imag, -1.0, places=15)
        self.assertAlmostEqual(res["quadrature_magnitude"], math.sqrt(2.0), places=15)
        self.assertAlmostEqual(res["quadrature_phase_deg"], -45.0, places=14)
        self.assertLess(res["destructive_magnitude"], 1e-15)
        self.assertAlmostEqual(res["recovered_group_delay_us"], 125.0, places=13)
        self.assertAlmostEqual(res["coherent_sum_db"], 20.0 * math.log10(4.0), places=14)
        self.assertAlmostEqual(res["array_white_noise_gain_db"], 10.0 * math.log10(4.0), places=14)
        self.assertGreater(res["pooled_energy_db"], res["block_average_of_db"] + 7.0)

    def test_e04_dft_parseval_and_zero_padded_linear_conv(self):
        res = cm.dft_parseval_and_zero_padded_conv_example()
        expected_dft = np.array([2.0 + 0.0j, 2.0 - 2.0j, -2.0 + 0.0j, 2.0 + 2.0j])
        np.testing.assert_allclose(res["dft_x"], expected_dft, rtol=0.0, atol=1e-14)
        self.assertAlmostEqual(res["time_energy"], 6.0, places=15)
        self.assertAlmostEqual(res["parseval_freq_energy"], 6.0, places=15)
        expected_lin = np.array([1.0, 2.5, -0.5, -1.5, 0.5, 0.0])
        np.testing.assert_allclose(res["linear_convolution"], expected_lin, rtol=0.0, atol=1e-15)
        np.testing.assert_allclose(res["circular_convolution_no_pad"], [1.5, 2.5, -0.5, -1.5], rtol=0.0, atol=1e-14)
        self.assertAlmostEqual(res["circular_aliasing_l2_error"], 0.5, places=14)
        self.assertLess(res["zero_padded_max_abs_error"], 1e-15)

    def test_e05_cola_and_wola_reconstruction(self):
        res = cm.stft_cola_window_reconstruction_example()
        self.assertLess(res["hann_ola_50pct"]["steady_max_abs_error"], 1e-14)
        self.assertLess(res["sqrt_hann_wola_50pct"]["steady_max_abs_error"], 1e-14)
        bad = res["unnormalized_squared_hann_50pct_counterexample"]
        self.assertAlmostEqual(bad["steady_cola_min"], 0.5, places=14)
        self.assertAlmostEqual(bad["steady_cola_max"], 1.0, places=14)
        self.assertAlmostEqual(bad["steady_cola_ripple"], 0.5, places=14)
        self.assertAlmostEqual(res["squared_hann_75pct"]["constant_cola_gain"], 1.5, places=14)
        self.assertLess(res["squared_hann_75pct"]["normalized_max_abs_error"], 1e-14)

    def test_e06_orthogonal_projection_and_four_subspaces(self):
        res = cm.orthogonal_projection_and_four_subspaces_example()
        expected_b_par = np.array([2.0 - 1.0j, -1.0 + 1.0j, -1.0 + 2.0j]) / 3.0
        expected_b_perp = np.array([1.0 + 1.0j, 1.0 - 1.0j, 1.0 + 1.0j]) / 3.0
        np.testing.assert_allclose(res["b_parallel"], expected_b_par, rtol=0.0, atol=1e-15)
        np.testing.assert_allclose(res["b_orthogonal"], expected_b_perp, rtol=0.0, atol=1e-15)
        self.assertAlmostEqual(res["gram_det"], 3.0, places=14)
        self.assertAlmostEqual(res["norm_sq_parallel"], 4.0 / 3.0, places=14)
        self.assertAlmostEqual(res["norm_sq_orthogonal"], 2.0 / 3.0, places=14)
        self.assertLess(res["idempotency_fro_error"], 1e-15)
        self.assertLess(res["orthogonality_inner_abs"], 1e-15)

    def test_e07_svd_condition_and_tikhonov_loading(self):
        res = cm.svd_condition_and_tikhonov_loading_example()
        self.assertAlmostEqual(res["condition_number_2"], 1000.0, places=12)
        np.testing.assert_allclose(res["unregularized_solution"], [1.0, 10.0], rtol=0.0, atol=1e-11)
        self.assertAlmostEqual(res["bound_tightness_ratio"], 1.0, places=11)
        expected_reg = np.array([100.0 / 100.1, 0.001 / 0.1001])
        np.testing.assert_allclose(res["tikhonov_solution"], expected_reg, rtol=0.0, atol=1e-12)
        self.assertGreater(res["error_reduction_factor"], 995.0)

    def test_e08_woodbury_and_schur_complement(self):
        res = cm.woodbury_and_schur_complement_example()
        self.assertAlmostEqual(res["sherman_morrison_denominator"], 7.0, places=15)
        self.assertLess(res["woodbury_vs_direct_fro_error"], 1e-15)
        self.assertAlmostEqual(res["interferer_attenuation_ratio"], 1.0 / 7.0, places=15)
        self.assertAlmostEqual(res["partition_det_a11"], 5.0, places=14)
        self.assertAlmostEqual(res["schur_complement_s"], 1.4, places=14)
        self.assertAlmostEqual(res["det_via_schur"], 7.0, places=14)

    def test_e09_rayleigh_and_gev_joint_diagonalization(self):
        res = cm.rayleigh_quotient_and_gev_diagonalization_example()
        np.testing.assert_allclose(res["r_nn_eigenvalues"], [1.0, 3.0], rtol=0.0, atol=1e-14)
        np.testing.assert_allclose(res["gev_eigenvalues"], [0.0, 4.0 / 3.0], rtol=0.0, atol=1e-14)
        self.assertAlmostEqual(res["max_gev_sinr"], 4.0 / 3.0, places=14)
        self.assertAlmostEqual(res["unwhitened_matched_sinr"], 1.0, places=14)
        self.assertAlmostEqual(res["whitening_sinr_gain_db"], 10.0 * math.log10(4.0 / 3.0), places=13)
        self.assertLess(res["b_whiten_fro_error"], 1e-14)
        self.assertLess(res["a_diag_fro_error"], 1e-14)

    def test_e10_vec_kronecker_and_khatri_rao(self):
        res = cm.vec_kronecker_and_khatri_rao_example()
        self.assertLess(res["vec_axb_max_abs_error"], 1e-14)
        self.assertGreater(res["wrong_hermitian_transpose_error"], 30.0)
        self.assertAlmostEqual(res["kronecker_trace"], 42.0, places=14)
        np.testing.assert_allclose(res["kronecker_eigenvalues"], [4.0, 8.0, 10.0, 20.0], rtol=0.0, atol=1e-14)
        self.assertEqual(res["difference_coarray_unique_lags"], [-3, -2, -1, 0, 1, 2, 3])
        self.assertLess(res["khatri_rao_vec_max_abs_error"], 1e-15)
        self.assertGreater(res["hadamard_psd_min_eigenvalue"], 0.3)

    def test_e11_wirtinger_mvdr_and_lcmv_hand_closed_forms(self):
        res = cm.wirtinger_lagrange_mvdr_and_lcmv_example()
        self.assertLess(res["wirtinger_fd_max_abs_error"], 1e-9)
        self.assertLess(res["directional_derivative_fd_error"], 1e-9)
        expected_mvdr = np.array([7.0 + 2.0j, 5.0 + 0.0j, 7.0 - 2.0j]) / 19.0
        expected_lcmv = np.array([3.0 + 1.0j, 2.0 + 0.0j, 3.0 - 1.0j]) / 8.0
        np.testing.assert_allclose(res["w_mvdr"], expected_mvdr, rtol=0.0, atol=1e-14)
        np.testing.assert_allclose(res["w_lcmv"], expected_lcmv, rtol=0.0, atol=1e-14)
        self.assertAlmostEqual(res["mvdr_output_power"], 7.0 / 19.0, places=14)
        self.assertAlmostEqual(res["lcmv_output_power"], 3.0 / 8.0, places=14)
        self.assertAlmostEqual(res["lcmv_minus_mvdr_power"], 1.0 / 152.0, places=14)

    def test_e12_nnls_kkt_and_mm_monotonicity(self):
        res = cm.nnls_kkt_and_mm_reweighted_l1_example()
        np.testing.assert_allclose(res["unconstrained_ls"], [22.0 / 15.0, -14.0 / 15.0], rtol=0.0, atol=1e-14)
        self.assertAlmostEqual(res["naive_clipped_objective"], 49.0 / 90.0, places=14)
        nnls = res["nnls_kkt"]
        np.testing.assert_allclose(nnls["q_opt"], [0.72, 0.0], rtol=0.0, atol=1e-14)
        np.testing.assert_allclose(nnls["dual_multiplier_mu"], [0.0, 0.42], rtol=0.0, atol=1e-14)
        self.assertAlmostEqual(nnls["objective_value"], 0.196, places=14)
        self.assertLess(nnls["complementary_slackness_max_abs"], 1e-14)
        self.assertTrue(res["mm_reweighted_l1"]["monotonically_non_increasing"])
        self.assertTrue(res["mm_reweighted_l1"]["surrogate_always_upper_bounds"])

    def test_e13_proper_complex_gaussian_dft_bins_and_crlb(self):
        res = cm.proper_complex_gaussian_and_crlb_example()
        self.assertTrue(res["proper_circular_case"]["is_proper_second_order"])
        self.assertFalse(res["improper_elliptical_case"]["is_proper_second_order"])
        self.assertAlmostEqual(res["improper_elliptical_case"]["pseudo_variance_abs"], math.sqrt(0.52), places=14)

        dft_p = res["real_dft_8pt_propriety"]
        self.assertAlmostEqual(dft_p["dc_pseudo_variance"], 1.0, places=14)
        self.assertAlmostEqual(dft_p["nyquist_pseudo_variance"], 1.0, places=14)
        self.assertLess(dft_p["interior_max_abs_pseudo_variance"], 1e-14)

        crlb = res["narrowband_tdoa_crlb"]
        expected_unknown_us = 1e6 / (2.0 * math.pi * 1000.0 * math.sqrt(250.0))
        self.assertAlmostEqual(crlb["crlb_std_unknown_phase_us"], expected_unknown_us, places=12)
        self.assertAlmostEqual(crlb["unknown_vs_known_variance_ratio"], 2.0, places=14)

    def test_e14_wiener_orthogonality_and_kalman_joseph(self):
        res = cm.wiener_orthogonality_and_kalman_joseph_example()
        w = res["wiener"]
        np.testing.assert_allclose(w["wiener_weights"], [4.0 / 11.0, 4.0 / 11.0], rtol=0.0, atol=1e-14)
        self.assertAlmostEqual(w["mmse"], 6.0 / 11.0, places=14)
        self.assertLess(w["orthogonality_residual_norm"], 1e-15)
        self.assertLess(w["mvdr_postfilter_equivalence_error"], 1e-15)

        kf = res["kalman"]
        np.testing.assert_allclose(kf["posterior_state"], [12.9375, 2.1875], rtol=0.0, atol=1e-14)
        np.testing.assert_allclose(
            kf["posterior_cov_joseph"], [[1.875, 0.375], [0.375, 1.875]], rtol=0.0, atol=1e-14
        )
        self.assertGreater(kf["perturbed_gain_naive_asymmetry"], 0.35)

    def test_e15_hungarian_matches_exhaustive_permutations(self):
        res = cm.hungarian_bipartite_assignment_example()
        self.assertAlmostEqual(res["greedy_total_cost"], 13.0, places=15)
        self.assertAlmostEqual(res["hungarian_total_cost"], 8.0, places=15)
        self.assertLess(res["duality_gap"], 1e-15)
        self.assertGreaterEqual(res["min_reduced_cost"], -1e-15)

        # Independent brute-force check on a 5x5 deterministic non-trivial matrix
        rng = np.random.default_rng(20260015)
        c5 = rng.uniform(-5.0, 20.0, size=(5, 5))
        brute_min = min(sum(c5[r, perm[r]] for r in range(5)) for perm in itertools.permutations(range(5)))
        opt5 = cm.hungarian_linear_sum_assignment(c5)
        self.assertAlmostEqual(opt5["total_cost"], brute_min, places=12)
        self.assertLess(opt5["duality_gap"], 1e-12)

    def test_e16_special_functions_independent_quadrature_and_identities(self):
        res = cm.bessel_legendre_and_spherical_harmonics_example()
        # Independent midpoint integral check for J_0(2), J_1(2), J_2(2)
        theta_mid = (np.arange(2000) + 0.5) * (math.pi / 2000.0)
        for order_n, actual_j in enumerate(res["bessel_j0_j1_j2_at_2"]):
            integral_j = float(np.mean(np.cos(order_n * theta_mid - 2.0 * np.sin(theta_mid))))
            self.assertAlmostEqual(actual_j, integral_j, places=13)
        self.assertLess(res["bessel_three_term_recurrence_error"], 1e-15)
        self.assertLess(res["jacobi_anger_abs_error"], 1e-14)
        self.assertLess(res["wronskian_identity_abs_error"], 1e-15)
        self.assertLess(res["open_sphere_b0_abs_at_pi"], 1e-14)
        self.assertAlmostEqual(
            res["rigid_sphere_b0_abs_at_pi"], res["rigid_sphere_b0_closed_form_abs_at_pi"], places=14
        )
        # Regression guard against the factor-of-pi simplification bug: 4*pi/sqrt(pi^2+1) == 4/sqrt(1+1/pi^2)
        self.assertAlmostEqual(
            res["rigid_sphere_b0_abs_at_pi"],
            4.0 * math.pi / math.sqrt(math.pi ** 2 + 1.0),
            places=14,
        )
        self.assertLess(res["legendre_orthogonality_fro_error"], 1e-14)
        self.assertLess(res["spherical_harmonic_addition_max_error"], 1e-14)

    def test_expanded_ch00_primitives_and_negative_controls(self):
        # 1. M-element ULA Dirichlet kernel vs direct geometric sum (including grating-lobe limit psi=2*pi)
        for psi in (0.0, 0.35 * math.pi, 0.5 * math.pi, 2.0 * math.pi, -1.25):
            dres = cm.dirichlet_ula_array_factor(4, psi)
            self.assertLess(dres["abs_error"], 1e-14)
        self.assertAlmostEqual(cm.dirichlet_ula_array_factor(4, 2.0 * math.pi)["normalized_magnitude"], 1.0, places=14)
        self.assertLess(cm.dirichlet_ula_array_factor(4, 0.5 * math.pi)["normalized_magnitude"], 1e-14)

        # 2. Hermitian square root R^{1/2} and R^{-1/2}
        r_mat = np.array([[2.0, 1.0 - 0.5j], [1.0 + 0.5j, 3.0]], dtype=np.complex128)
        sq = cm.hermitian_matrix_sqrt(r_mat)
        self.assertLess(sq["square_fro_error"], 1e-14)
        self.assertLess(sq["whitening_fro_error"], 1e-14)

        # 3. Matrix Wirtinger log-determinant gradient: d ln|det W| / d W^* = 0.5 W^{-H}
        w_mat = np.array([[2.0 + 1.0j, 1.0], [0.5j, 1.0 - 1.0j]], dtype=np.complex128)
        lw = cm.matrix_wirtinger_logdet_check(w_mat)
        self.assertLess(lw["max_abs_error"], 1e-9)

        # 4. Gauss-Markov GLS via Cholesky whitening & BLUE optimality
        a_gls = np.array([[1.0], [1.0]], dtype=np.complex128)
        b_gls = np.array([1.2, 0.8], dtype=np.complex128)
        c_gls = np.array([[1.0, 0.8], [0.8, 2.0]], dtype=np.complex128)
        gls = cm.gauss_markov_gls_solve(a_gls, b_gls, c_gls)
        self.assertLess(gls["weighted_orthogonality_norm"], 1e-14)
        self.assertGreaterEqual(gls["cov_ols_minus_gls_min_eig"], 0.0)
        # Correlated strong noise produces negative weight on noisier mic when C = [[1, 1.2], [1.2, 2]]
        c_neg = np.array([[1.0, 1.2], [1.2, 2.0]], dtype=np.complex128)
        w_gls_neg = np.linalg.solve(c_neg, a_gls)
        w_gls_neg = (w_gls_neg / (a_gls.conj().T @ w_gls_neg)).ravel().real
        self.assertGreater(w_gls_neg[0], 1.0)
        self.assertLess(w_gls_neg[1], 0.0)

        # 5. Binomial zero-failure upper bound & Wilson score interval
        self.assertAlmostEqual(cm.binomial_zero_failure_upper_bound(299, 0.05), 1.0 - 0.05 ** (1.0 / 299.0), places=15)
        self.assertLess(cm.binomial_zero_failure_upper_bound(299, 0.05), 0.01)
        self.assertGreater(cm.binomial_zero_failure_upper_bound(298, 0.05), 0.01)
        w_ci = cm.wilson_score_interval(0, 20, z_value=1.959963984540054)
        self.assertEqual(w_ci["lower"], 0.0)
        self.assertGreater(w_ci["upper"], 0.16)

        # 6. First-order real spherical harmonics ACN / SN3D / N3D / orthonormal
        foa = cm.real_spherical_harmonic_foa(math.pi / 3.0, math.pi / 4.0)
        np.testing.assert_allclose(foa["acn_n3d"][1:], math.sqrt(3.0) * foa["acn_sn3d"][1:], rtol=0.0, atol=1e-15)
        np.testing.assert_allclose(foa["acn_orthonormal"], foa["acn_n3d"] / math.sqrt(4.0 * math.pi), rtol=0.0, atol=1e-15)
        # Addition theorem check on real orthonormal n=1 harmonics: sum_m Y^2 = 3/(4*pi)
        self.assertAlmostEqual(float(np.sum(foa["acn_orthonormal"][1:] ** 2)), 3.0 / (4.0 * math.pi), places=15)

        # Strict negative controls on new primitives
        with self.assertRaises(ValueError):
            cm.dirichlet_ula_array_factor(True, 0.5)
        with self.assertRaises(ValueError):
            cm.hermitian_matrix_sqrt(np.diag([1.0, -1.0]))
        with self.assertRaises(ValueError):
            cm.matrix_wirtinger_logdet_check(np.array([[1.0, 1.0], [1.0, 1.0]]))
        with self.assertRaises(ValueError):
            cm.wilson_score_interval(21, 20)
        with self.assertRaises(ValueError):
            cm.binomial_zero_failure_upper_bound(10, 1.0)

    def test_strict_validation_rejects_booleans_complex_nan_and_invalid_domains(self):
        # Boolean rejected where int or float is required
        with self.assertRaises(ValueError):
            cm.phasor_delay_and_decibel_budget_example(num_mics=True)
        with self.assertRaises(ValueError):
            cm.coordinate_convention_bridge(True, 0.0, 1.0)
        with self.assertRaises(ValueError):
            cm.dft_parseval_and_zero_padded_conv_example(x=[True, False, True, False])
        with self.assertRaises(ValueError):
            cm.complex_column_space_projection(np.array([[True], [False]]), np.array([1.0, 0.0]))

        # Complex rejected for real-valued APIs without discarding imaginary parts
        with self.assertRaises(ValueError):
            cm.so3_rotation_matrix_composition_example(vector=[1.0 + 0.0j, 2.0, 2.0])
        with self.assertRaises(ValueError):
            cm.hungarian_linear_sum_assignment(np.array([[1.0 + 1.0j, 2.0], [3.0, 4.0]]))

        # NaN / Inf rejected
        with self.assertRaises(ValueError):
            cm.narrowband_tdoa_crlb(float("nan"), 10.0, 5)
        with self.assertRaises(ValueError):
            cm.tikhonov_svd_solve(np.eye(2), np.array([1.0, float("inf")]), 0.1)

        # Non-Hermitian and indefinite matrices rejected in Cholesky/LCMV
        with self.assertRaises(ValueError):
            cm.solve_lcmv_weights(np.array([[1.0, 2.0], [0.0, 1.0]]), np.ones((2, 1)), np.ones(1))
        with self.assertRaises(ValueError):
            cm.solve_hermitian_gev(np.eye(2), np.diag([1.0, -0.5]))

        # Domain checks
        with self.assertRaises(ValueError):
            cm.legendre_polynomial(2, 1.05)
        with self.assertRaises(ValueError):
            cm.spherical_bessel_jy(0, 0.0)
        with self.assertRaises(ValueError):
            cm.fresnel_spherical_wave_approx(0.10, 0.10, 0.0, 1000.0, 0.20)

    def test_scale_invariant_pd_large_argument_bessel_and_singular_guards(self):
        # 1. Scale-invariant Hermitian / strictly PD checks across 24 orders of magnitude
        gev_tiny = cm.solve_hermitian_gev(2e-12 * np.eye(2), 1e-12 * np.eye(2))
        np.testing.assert_allclose(gev_tiny["gev_eigenvalues"], [2.0, 2.0], rtol=0.0, atol=1e-13)
        w_tiny = cm.solve_lcmv_weights(1e-12 * np.eye(2), np.array([[1.0], [0.0]]), np.array([1.0]))
        np.testing.assert_allclose(w_tiny, [1.0 + 0.0j, 0.0 + 0.0j], rtol=0.0, atol=1e-14)
        w_huge = cm.solve_lcmv_weights(1e12 * np.eye(2), np.array([[1.0], [0.0]]), np.array([1.0]))
        np.testing.assert_allclose(w_huge, [1.0 + 0.0j, 0.0 + 0.0j], rtol=0.0, atol=1e-14)
        with self.assertRaises(ValueError):
            cm.solve_lcmv_weights(np.zeros((2, 2)), np.array([[1.0], [0.0]]), np.array([1.0]))
        with self.assertRaises(ValueError):
            cm.solve_lcmv_weights(
                np.array([[1e-12, 2e-12], [0.0, 1e-12]]),
                np.array([[1.0], [0.0]]),
                np.array([1.0]),
            )

        # 2. Large-argument Bessel quadrature & upward recurrence free of alternating-series cancellation
        j0_30 = cm.cylindrical_bessel_j(0, 30.0)
        j1_30 = cm.cylindrical_bessel_j(1, 30.0)
        j2_30 = cm.cylindrical_bessel_j(2, 30.0)
        self.assertLess(abs(j0_30 + j2_30 - (2.0 / 30.0) * j1_30), 1e-14)
        self.assertAlmostEqual(cm.cylindrical_bessel_j(-1, 25.0), -cm.cylindrical_bessel_j(1, 25.0), places=15)

        x_lg = 25.0
        sj5, sy5 = cm.spherical_bessel_jy(5, x_lg)
        dj5, dy5 = cm.spherical_bessel_derivative(5, x_lg)
        wronskian_25 = sj5 * dy5 - dj5 * sy5
        self.assertLess(abs(wronskian_25 - 1.0 / (x_lg * x_lg)), 1e-14)
        sj_mix, sy_mix = cm.spherical_bessel_jy(4, np.array([0.5, 25.0]))
        self.assertAlmostEqual(float(sj_mix[0]), cm.spherical_bessel_jy(4, 0.5)[0], places=15)
        self.assertAlmostEqual(float(sj_mix[1]), cm.spherical_bessel_jy(4, 25.0)[0], places=15)
        self.assertAlmostEqual(float(sy_mix[1]), cm.spherical_bessel_jy(4, 25.0)[1], places=15)

        # 3. Underdetermined / singular input guards in tikhonov_svd_solve & woodbury_inverse_update
        with self.assertRaises(ValueError):
            cm.tikhonov_svd_solve(np.array([[1.0, 2.0]]), np.array([1.0]), 0.0)
        with self.assertRaises(ValueError):
            cm.tikhonov_svd_solve(np.array([[1.0, 1.0], [1.0, 1.0]]), np.array([1.0, 1.0]), 0.0)
        with self.assertRaises(ValueError):
            cm.woodbury_inverse_update(np.zeros((2, 2)), np.ones((2, 1)), np.eye(1), np.ones((1, 2)))
        with self.assertRaises(ValueError):
            cm.woodbury_inverse_update(np.eye(2), np.ones((2, 1)), np.zeros((1, 1)), np.ones((1, 2)))
        with self.assertRaises(ValueError):
            cm.woodbury_inverse_update(
                np.eye(2),
                np.array([[1.0], [0.0]]),
                np.array([[-1.0]]),
                np.array([[1.0, 0.0]]),
            )

    def test_regression_small_magnitude_bessel_length1_dft_and_underdetermined_tikhonov(self):
        # 1. Elementwise relative convergence for small-magnitude Bessel series: J_35(6.0) & j_10(0.1)
        with decimal.localcontext() as ctx:
            ctx.prec = 60
            half_x = decimal.Decimal("3")
            ref_j35 = float(
                sum(
                    ((-1) ** k) * (half_x ** (2 * k + 35))
                    / decimal.Decimal(math.factorial(k) * math.factorial(35 + k))
                    for k in range(60)
                )
            )
            dx = decimal.Decimal("0.1")
            ref_sj10 = float(
                sum(
                    ((-1) ** k) * (dx ** (2 * k + 10))
                    / decimal.Decimal(
                        (2 ** k)
                        * math.factorial(k)
                        * (math.factorial(20 + 2 * k + 1) // ((2 ** (10 + k)) * math.factorial(10 + k)))
                    )
                    for k in range(40)
                )
            )

        got_j35 = cm.cylindrical_bessel_j(35, 6.0)
        self.assertLess(abs(got_j35 - ref_j35) / abs(ref_j35), 1e-14)
        got_j35_arr = cm.cylindrical_bessel_j(35, np.array([0.0, 6.0]))
        self.assertEqual(float(got_j35_arr[0]), 0.0)
        self.assertLess(abs(float(got_j35_arr[1]) - ref_j35) / abs(ref_j35), 1e-14)
        got_j0_zero = cm.cylindrical_bessel_j(0, np.array([0.0, 6.0]))
        self.assertAlmostEqual(float(got_j0_zero[0]), 1.0, places=15)

        got_sj10, _ = cm.spherical_bessel_jy(10, 0.1)
        self.assertLess(abs(got_sj10 - ref_sj10) / abs(ref_sj10), 1e-14)
        # Mixed small/large array exercises the second small_mask power-series branch
        got_sj10_mix, _ = cm.spherical_bessel_jy(10, np.array([0.1, 15.0]))
        self.assertLess(abs(float(got_sj10_mix[0]) - ref_sj10) / abs(ref_sj10), 1e-14)
        self.assertAlmostEqual(float(got_sj10_mix[1]), cm.spherical_bessel_jy(10, 15.0)[0], places=15)

        # 2. Length-1 DFT Parseval & zero-padded convolution guard (empty interior k_int)
        res_len1 = cm.dft_parseval_and_zero_padded_conv_example(x=[3.0], h=[2.0], n_pad=1)
        self.assertAlmostEqual(res_len1["time_energy"], 9.0, places=15)
        self.assertAlmostEqual(res_len1["parseval_freq_energy"], 9.0, places=15)
        self.assertEqual(res_len1["conjugate_symmetry_error"], 0.0)
        np.testing.assert_allclose(res_len1["linear_convolution"], [6.0], rtol=0.0, atol=1e-15)
        self.assertLess(res_len1["zero_padded_max_abs_error"], 1e-15)

        res_len1_h2 = cm.dft_parseval_and_zero_padded_conv_example(x=[3.0], h=[2.0, -1.0], n_pad=2)
        self.assertEqual(res_len1_h2["conjugate_symmetry_error"], 0.0)
        np.testing.assert_allclose(res_len1_h2["linear_convolution"], [6.0, -3.0], rtol=0.0, atol=1e-15)

        # 3. Underdetermined tikhonov_svd_solve with tiny eps=1e-20 where normal equations are singular
        res_ud = cm.tikhonov_svd_solve(np.array([[1.0, 2.0]]), np.array([1.0]), loading_eps=1e-20)
        np.testing.assert_allclose(res_ud["solution"], [0.2 + 0.0j, 0.4 + 0.0j], rtol=0.0, atol=1e-14)
        self.assertTrue(math.isinf(res_ud["condition_number_2"]))
        self.assertTrue(math.isnan(res_ud["normal_eq_consistency_error"]))


class Chapter00ExercisesCLITests(unittest.TestCase):
    def test_in_process_registry_and_strict_json(self):
        self.assertEqual(len(ex00.EXERCISE_IDS), 16)
        payload = ex00.run_experiments()
        encoded = json.dumps(payload, ensure_ascii=False, allow_nan=False)
        parsed = strict_json_loads(encoded)
        self.assertEqual(set(parsed["exercises"]), set(ex00.EXERCISE_IDS))

    def test_direct_script_cli_modes(self):
        script = ROOT / "codes/chapters/ch00/chapter00_exercises.py"
        proc_list = subprocess.run(
            [sys.executable, str(script), "--list"],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
            check=True,
        )
        lines = [ln.strip() for ln in proc_list.stdout.splitlines() if ln.strip()]
        self.assertEqual(len(lines), 16)
        self.assertTrue(lines[0].startswith("E00-01:"))
        self.assertTrue(lines[-1].startswith("E00-16:"))

        proc_list_json = subprocess.run(
            [sys.executable, str(script), "--list", "--json"],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
            check=True,
        )
        list_payload = strict_json_loads(proc_list_json.stdout)
        self.assertEqual(list_payload["exercise_ids"], list(ex00.EXERCISE_IDS))

        proc_single = subprocess.run(
            [sys.executable, str(script), "--exercise", "E00-11", "--json"],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
            check=True,
        )
        single_payload = strict_json_loads(proc_single.stdout)
        self.assertEqual(list(single_payload["exercises"]), ["E00-11"])
        self.assertAlmostEqual(single_payload["exercises"]["E00-11"]["mvdr_output_power"], 7.0 / 19.0, places=14)


if __name__ == "__main__":
    unittest.main()
