"""Independent algebra, input boundaries and update-time contracts for ch13."""
import json
import math
import unittest
import unittest.mock
import numpy as np
from codes.chapters.ch13.core.distributed import (
    known_models, mwf_weights, mse_components, compressed_mwf, distributed_updates,
    jacobi_control, step_size_control, tree_sum_control, gevd_control,
)
from codes.chapters.ch13.core.dasf_and_ptp import (
    ptp_vs_sro_phase_budget, centralized_lcmv, dasf_lcmv_updates, dasf_unification_control,
)
from codes.chapters.ch13.chapter13_exercises import (
    run_experiments, resource_control, dwacd_sro_gating_control,
)


class DistributedTests(unittest.TestCase):
    def setUp(self):
        self.a = np.array([1., .5, 2., -.5])
        self.rs = np.outer(self.a, self.a)
        self.rn = np.eye(4)
        self.rn[0, 2] = self.rn[2, 0] = .2
        self.rn[0, 3] = self.rn[3, 0] = .8

    def test_white_and_correlated_rational_anchors(self):
        for noise, expected, mse in [(np.eye(4), np.array([2, 1, 4, -1])/13, 2/13),
                                     (self.rn, np.array([25, 4, 11, -24])/69, 8/69)]:
            w = mwf_weights(self.rs, noise)
            np.testing.assert_allclose(w, expected, atol=5e-15)
            self.assertAlmostEqual(mse_components(self.rs, noise, w)['total_mse'], mse)
            # Independent expansion for a nonoptimal half-sized weight.
            q = expected/2
            scalar = (q@self.a-1)**2+q@noise@q
            self.assertAlmostEqual(mse_components(self.rs, noise, q)['total_mse'], scalar)

    def test_complete_arbitrary_weight_cost(self):
        result = mse_components(np.diag([2., 1.]), np.eye(2), [2., 3.])
        self.assertEqual(result['target_distortion'], 11.)
        self.assertEqual(result['noise_power'], 13.)
        self.assertEqual(result['total_mse'], 24.)
        zero = mse_components(self.rs, np.eye(4), np.zeros(4))
        self.assertEqual(zero['total_mse'], 1.)

    def test_complex_reference_conjugation(self):
        a = np.array([1., .5j, 2., -.5])
        rs = np.outer(a, a.conj())
        expected = np.array([-1j, .5, -2j, .5j])/13
        w = mwf_weights(rs, np.eye(4), 1)
        np.testing.assert_allclose(w, expected, atol=2e-16)
        self.assertAlmostEqual(np.vdot(w, a).imag, 11/26)
        self.assertAlmostEqual(mse_components(rs, np.eye(4), w, 1)['total_mse'], 1/26)

    def test_compression_and_independent_kkt(self):
        t = np.array([[1., 0., 0., 0.], [0., 1., 0., 0.], [0., 0., 2., -.5]])
        result = compressed_mwf(self.rs, self.rn, t)
        normal = np.array([0., 0., .5, 2.])
        kkt = np.block([[self.rs+self.rn, normal[:, None]], [normal[None, :], np.zeros((1, 1))]])
        independent = np.linalg.solve(kkt, np.r_[self.a, 0])[:4]
        np.testing.assert_allclose(result['weights'], independent, atol=5e-16)
        self.assertAlmostEqual(result['components']['total_mse'], 2/13)
        true_w = np.array([25, 4, 11, -24])/69
        t[2, 2:] = true_w[2:]
        np.testing.assert_allclose(compressed_mwf(self.rs, self.rn, t)['weights'], true_w, atol=5e-16)

    def test_normalized_nonzero_scale_and_complex_phase(self):
        expected = np.array([2, 1, 4, -1])/13
        for scale in (1e-300, 1e300, 2j):
            t = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 2, -.5]], complex)
            t[2] *= scale
            result = compressed_mwf(self.rs, np.eye(4), t)
            np.testing.assert_allclose(result['weights'], expected, atol=3e-15)
            np.testing.assert_allclose(result['normalized_projection'].conj().T@result['compressed_weights'], expected, atol=3e-15)

    def test_complex_projection_uses_conjugated_weights(self):
        a = np.array([1., .5j, 2j, -.5])
        rs = np.outer(a, a.conj())
        t = np.zeros((3, 4), complex); t[:2, :2] = np.eye(2); t[2, 2:] = a[2:].conj()
        expected = a/6.5
        np.testing.assert_allclose(compressed_mwf(rs, np.eye(4), t)['weights'], expected, atol=5e-16)

    def test_invalid_projection_shapes_zero_and_dependent_rows(self):
        for t in (np.zeros((3, 4)), np.ones((2, 4)), np.eye(5), np.ones(4), [[True]*4]):
            with self.assertRaises((ValueError, np.linalg.LinAlgError)):
                compressed_mwf(self.rs, np.eye(4), t)

    def test_covariance_validation_precedes_solve(self):
        for rs, rn in [(np.eye(2), [[1, 2], [2, 1]]),
                       ([[1, 1], [0, 1]], np.eye(2)), (np.zeros((2, 2)), np.zeros((2, 2))),
                       (np.diag([1., 0.]), np.diag([1., 1e-14])), (np.eye(2), np.eye(3))]:
            with self.assertRaises((ValueError, np.linalg.LinAlgError)):
                mwf_weights(rs, rn)
        for bad in (True, -1, 4, 1.5):
            with self.assertRaises(ValueError):
                mwf_weights(self.rs, np.eye(4), bad)

    def test_covariance_scale_invariance_and_physical_units(self):
        expected = np.array([25, 4, 11, -24])/69
        for scale in (1e-300, 1e300):
            w = mwf_weights(self.rs*scale, self.rn*scale)
            np.testing.assert_allclose(w, expected, atol=5e-15)
            result = mse_components(self.rs*scale, self.rn*scale, w)
            self.assertAlmostEqual(result['total_mse']/scale, 8/69)
        with self.assertRaises(ValueError):
            mse_components(np.eye(2), np.eye(2), [1e308, 1e308])

    def test_zero_reference_is_explicit(self):
        result = mse_components(np.diag([1., 0.]), np.eye(2), [0., 0.], 1)
        self.assertEqual(result['reference_target_power'], 0.)
        self.assertIsNone(result['normalized_mse'])

    def _updates(self, **kwargs):
        return distributed_updates(self.rs, self.rn, ((0, 1), (2, 3)), (0, 2),
                                   [self.a[:2], self.a[2:]], **kwargs)

    def test_round_robin_preserves_stale_receivers_and_rebuilds_current_outputs(self):
        result = self._updates(max_updates=6, tolerance=0.)
        history = result['history']
        expected_j = [2/13, .1251372118551044, .1180174009356284,
                      .11605609650078663, .11594737978282565, .1159422132581912]
        for i, entry in enumerate(history):
            k = i % 2
            self.assertEqual(entry['active_nodes'], [k])
            self.assertEqual(entry['total_solve_count'], i+3)
            self.assertAlmostEqual(entry['solutions'][0]['components']['normalized_mse'], expected_j[i], places=12)
            if i:
                np.testing.assert_array_equal(entry['receiver_coefficients_raw_coordinates'][1-k],
                                              history[i-1]['receiver_coefficients_raw_coordinates'][1-k])
                np.testing.assert_array_equal(entry['outputs_at_last_solve'][1-k], history[i-1]['outputs_at_last_solve'][1-k])
                self.assertGreater(entry['output_age_in_update_solves'][1-k], 0)
            for peer in (0, 1):
                expected = entry['current_receiver_projections'][peer].conj().T@entry['receiver_coefficients_raw_coordinates'][peer]
                np.testing.assert_array_equal(entry['cached_outputs'][peer], expected)
        saved = history[0]['compressions_after'][0].copy()
        result['compressions'][0][:] = 100
        np.testing.assert_array_equal(history[0]['compressions_after'][0], saved)

    def test_independent_white_two_step_current_broadcast_rational_anchor(self):
        result = distributed_updates(self.rs, np.eye(4), ((0, 1), (2, 3)), (0, 2),
                                     [[1., 0.], [1., 0.]], max_updates=2, tolerance=0.)
        first, second = result['history']
        np.testing.assert_allclose(first['cached_outputs'][0], np.array([4, 2, 8, 0])/25, atol=1e-16)
        np.testing.assert_allclose(second['cached_outputs'][0], np.array([52, 26, 64, -16])/325, atol=2e-16)
        np.testing.assert_array_equal(first['receiver_coefficients_raw_coordinates'][0], second['receiver_coefficients_raw_coordinates'][0])
        np.testing.assert_array_equal(first['outputs_at_last_solve'][0], second['outputs_at_last_solve'][0])
        self.assertAlmostEqual(second['cached_components'][0]['total_mse'], 23108/105625, places=14)
        self.assertEqual(second['total_solve_count'], 4)

    def test_simultaneous_snapshot_and_budget(self):
        result = self._updates(schedule='simultaneous', max_updates=3, tolerance=0.)
        self.assertEqual(result['update_solve_count'], 2)
        self.assertEqual(result['unused_update_budget'], 1)
        self.assertEqual(result['total_solve_count'], 4)
        entry = result['history'][0]
        self.assertEqual(entry['active_nodes'], [0, 1])
        for k in (0, 1):
            for peer in (0, 1):
                np.testing.assert_array_equal(entry['solutions'][k]['incoming_compressions'][peer], entry['compressions_before'][peer])
            expected = entry['current_receiver_projections'][k].conj().T@entry['receiver_coefficients_raw_coordinates'][k]
            np.testing.assert_array_equal(entry['cached_outputs'][k], expected)
            np.testing.assert_array_equal(entry['solutions'][k]['effective_weights_after_broadcast'], expected)
        self.assertFalse(np.allclose(entry['cached_outputs'][0], entry['solutions'][0]['weights']))

    def test_actual_residual_stopping_and_count(self):
        result = self._updates(max_updates=40, tolerance=1e-8)
        self.assertEqual(result['status'], 'known_covariance_residual_reached')
        self.assertEqual(result['update_solve_count'], 14)
        self.assertEqual(result['total_solve_count'], 16)
        self.assertLessEqual(max(result['history'][-1]['normal_equation_relative_residuals']), 1e-8)
        self.assertGreater(max(result['history'][-2]['normal_equation_relative_residuals']), 1e-8)

    def test_zero_peer_only_explicit_local_policy(self):
        arguments = (self.rs, np.eye(4), ((0, 1), (2, 3)), (0, 2), [np.zeros(2), np.zeros(2)])
        with self.assertRaises(np.linalg.LinAlgError):
            distributed_updates(*arguments)
        result = distributed_updates(*arguments, zero_policy='local', max_updates=1, tolerance=0.)
        solve = result['history'][0]['solutions'][0]
        self.assertEqual(solve['missing_peers'], [1])
        np.testing.assert_allclose(solve['weights'], [4/9, 2/9, 0, 0])
        self.assertEqual(result['receiver_peer_maps'][1], [])
        np.testing.assert_allclose(result['cached_outputs'][1], [0, 0, 16/21, -4/21])
        for kwargs in ({'relaxation': 0}, {'max_updates': True}, {'tolerance': -1}, {'schedule': 'unknown'}):
            with self.assertRaises(ValueError):
                self._updates(**kwargs)

    def test_target_rank_counterexample(self):
        a = np.array([[1, 0], [0, 1], [1, 1], [1, -1]], float)
        rs = a@a.T; rn = np.eye(4); t = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 1]]
        full = np.column_stack([mwf_weights(rs, rn, r) for r in (0, 1)])
        np.testing.assert_allclose(full, a/4)
        self.assertEqual(np.linalg.matrix_rank(full[2:]), 2)
        self.assertAlmostEqual(compressed_mwf(rs, rn, t, 0)['components']['total_mse'], 1/4)
        self.assertAlmostEqual(compressed_mwf(rs, rn, t, 1)['components']['total_mse'], 1/2)

    def test_jacobi_controls(self):
        bad = jacobi_control(); good = jacobi_control(.5, 40)
        np.testing.assert_allclose(bad['trace'][1:, 0], [1., -.8, 2.44, -3.392], atol=1e-14)
        self.assertAlmostEqual(bad['spectral_radius'], 1.8)
        self.assertAlmostEqual(good['spectral_radius'], .95)
        # Symmetric rhs excites only the all-ones eigenvector, multiplier -.4.
        np.testing.assert_allclose(good['trace'][-1], np.full(3, 5/14), atol=1e-14)

    def test_infinite_step_conditions_are_not_empirical(self):
        result = step_size_control(8)
        np.testing.assert_allclose(result['harmonic']['finite_values'], [1, 1/2, 1/3, 1/4, 1/5, 1/6, 1/7, 1/8])
        self.assertTrue(result['harmonic']['limit_zero'])
        self.assertFalse(result['constant_half']['limit_zero'])
        self.assertEqual(result['geometric_half']['infinite_sum'], '1')
        self.assertEqual(result['geometric_half']['infinite_square_sum'], '1/3')

    def test_actual_two_pass_tree_messages(self):
        result = tree_sum_control()
        self.assertEqual([(m['from'], m['to'], m['value']) for m in result['edge_messages']],
                         [(3, 2, 3), (2, 1, 5), (1, 2, 6), (2, 3, 6)])
        np.testing.assert_array_equal(result['totals'], [6, 6, 6])
        np.testing.assert_array_equal(result['subtract_own'], [5, 4, 3])
        np.testing.assert_array_equal(result['neighbors_self_only'], [2, 4, 2])

    def test_gevd_independent_analytic_values_and_metric(self):
        result = gevd_control()
        np.testing.assert_allclose(result['generalized_eigenvalues'], [5, 1], atol=3e-15)
        np.testing.assert_allclose(result['rank1_weights'], [0, 2/5], atol=5e-16)
        self.assertAlmostEqual(result['components']['total_mse'], 1/5)
        np.testing.assert_allclose(result['noise_metric_gram'], np.eye(2), atol=8e-16)
        np.testing.assert_allclose(result['generalized_residual'], 0, atol=3e-15)
        np.testing.assert_allclose(result['ordinary_eigenvalues'], [(7+np.sqrt(34))/2, (7-np.sqrt(34))/2])

    def test_resource_payload_arithmetic(self):
        result = resource_control()
        self.assertEqual(result['two_nodes_PCM16_raw_upload_bps'], 1024000)
        self.assertEqual(result['two_nodes_PCM16_scalar_exchange_bps'], 512000)
        self.assertEqual(result['two_nodes_float32_scalar_exchange_bps'], 1024000)
        self.assertEqual(result['two_nodes_STFT_complex_float32_exchange_bps'], 4112000)
        self.assertEqual(result['network'][1]['compressed_full_mesh_unicast_PCM16_bits_per_second'], 3072000)
        self.assertEqual(result['network'][1]['network_full_complex128_covariance_bytes'], 1600)

    def test_ptp_vs_sro_phase_budget_hand_anchors(self):
        budget = ptp_vs_sro_phase_budget()
        self.assertAlmostEqual(budget['sample_period_us'], 125/6, places=13)
        ptp = budget['ptp_hardware_sync']
        self.assertAlmostEqual(ptp['single_node_samples'], .048, places=15)
        self.assertAlmostEqual(ptp['single_node_phase_deg'], .864, places=14)
        self.assertAlmostEqual(ptp['worst_two_node_samples'], .096, places=15)
        self.assertAlmostEqual(ptp['worst_two_node_phase_deg'], 1.728, places=14)
        self.assertAlmostEqual(ptp['equal_sum_power_gain'], math.cos(.0048*math.pi)**2, places=15)
        self.assertAlmostEqual(ptp['target_waveform_phase_mse'], math.sin(.0048*math.pi)**2, places=18)
        self.assertFalse(ptp['accumulates_over_time'])

        sro = budget['free_running_sro']
        self.assertAlmostEqual(sro['sample_drift_per_s'], 4.8, places=14)
        self.assertAlmostEqual(sro['seconds_per_integer_sample_slip'], 5/24, places=15)
        self.assertAlmostEqual(sro['phase_drift_deg_per_s'], 86.4, places=13)
        self.assertAlmostEqual(sro['first_destructive_null_s'], 25/12, places=14)
        self.assertAlmostEqual(sro['first_destructive_null_samples'], 10., places=14)
        self.assertAlmostEqual(sro['wola']['drift_samples_per_hop'], .0512, places=15)
        self.assertAlmostEqual(sro['wola']['phase_drift_deg_per_hop'], .9216, places=14)
        self.assertAlmostEqual(sro['wola']['hops_per_integer_sample_slip'], 625/32, places=13)
        with self.assertRaises(ValueError):
            ptp_vs_sro_phase_budget(fs_hz=True)
        # Quarter-period PTP bound yields destructive coherent-sum zero/near-zero without domain error.
        with unittest.mock.patch("codes.chapters.ch13.core.dasf_and_ptp.math.cos", return_value=0.0):
            null_budget = ptp_vs_sro_phase_budget()
        self.assertEqual(null_budget['ptp_hardware_sync']['equal_sum_attenuation_db'], math.inf)

    def test_dasf_unification_stale_subspace_and_dof_lockup(self):
        ctrl = dasf_unification_control()
        r1 = ctrl['rank1_shared_subspace']
        np.testing.assert_allclose(r1['r_nn_inv_a'], np.array([25, 4, 11, -24])/8, atol=2e-15)
        self.assertAlmostEqual(r1['q_snr'], 61/8, places=14)
        np.testing.assert_allclose(r1['mwf_central_weights'], np.array([25, 4, 11, -24])/69, atol=2e-15)
        np.testing.assert_allclose(r1['mvdr_central_weights'], np.array([25, 4, 11, -24])/61, atol=2e-15)
        self.assertAlmostEqual(r1['mvdr_min_noise_power'], 8/61, places=15)
        self.assertAlmostEqual(r1['gevd_lambda_max'], 69/8, places=14)
        np.testing.assert_allclose(
            r1['gevd_central_weights'],
            np.array([25, 4, 11, -24])/(2*math.sqrt(122)),
            atol=2e-15,
        )
        np.testing.assert_allclose(r1['compressed_node1_mvdr_h'], np.array([25, 4, 1])/61, atol=2e-15)
        self.assertLess(r1['compressed_node1_mwf_error_l2'], 1e-14)
        self.assertLess(r1['compressed_node1_mvdr_error_l2'], 1e-14)
        self.assertLess(r1['compressed_node1_gevd_error_l2'], 1e-14)

        lcmv2 = ctrl['two_constraint_lcmv']
        expected_w2 = np.array([125, -13, 13, -117])/203
        np.testing.assert_allclose(lcmv2['central_weights'], expected_w2, atol=2e-15)
        self.assertAlmostEqual(lcmv2['central_noise_power'], 34/203, places=15)
        np.testing.assert_allclose(lcmv2['stale_mvdr_subspace_h'], np.array([5210, -2002, 182])/10397, atol=2e-15)
        np.testing.assert_allclose(lcmv2['stale_mvdr_subspace_weights'], np.array([5210, -2002, 2002, -4368])/10397, atol=2e-15)
        self.assertAlmostEqual(lcmv2['stale_mvdr_subspace_noise_power'], 2116/10397, places=15)
        self.assertAlmostEqual(lcmv2['stale_excess_noise_power'], 76050/2110591, places=15)
        self.assertTrue(lcmv2['dasf_r1_monotone'])
        np.testing.assert_allclose(lcmv2['dasf_r1_final_weights'], expected_w2, atol=1e-12)
        self.assertAlmostEqual(lcmv2['dasf_r1_final_noise_power'], 34/203, places=13)

        dof3 = ctrl['three_constraint_dof_boundary']
        expected_w3 = np.array([34, -26, 26, 0])/73
        np.testing.assert_allclose(dof3['central_weights'], expected_w3, atol=2e-15)
        self.assertAlmostEqual(dof3['central_noise_power'], 196/365, places=14)
        self.assertEqual(dof3['r1_locked_remaining_dofs'], [0, 0])
        np.testing.assert_allclose(dof3['r1_locked_weights'], [1., 0., 0., 0.], atol=2e-15)
        self.assertAlmostEqual(dof3['r1_locked_noise_power'], 1., places=15)
        self.assertEqual(dof3['r2_unlocked_remaining_dofs'], [1, 1])
        np.testing.assert_allclose(dof3['r2_unlocked_weights'], expected_w3, atol=1e-13)
        self.assertAlmostEqual(dof3['r2_unlocked_noise_power'], 196/365, places=14)

        for bad_call in (
            lambda: centralized_lcmv([[1., 2.], [2., 1.]], [[1.], [0.]], [1.]),
            lambda: centralized_lcmv(np.eye(2), [[1., 1.], [1., 1.]], [1., 0.]),
            lambda: centralized_lcmv(np.eye(2), [[1.], [0.]], [1.], condition_limit=0.0),
            lambda: centralized_lcmv(np.eye(2), [[1.], [0.]], [1.], condition_limit=True),
            lambda: centralized_lcmv(np.eye(2), [[1.], [0.]], [1.], condition_limit=float('nan')),
            lambda: dasf_lcmv_updates(self.rn, np.eye(4), [1., 0., 0., 0.], node_sizes=(2, 2), num_filters=1),
            lambda: dasf_lcmv_updates(self.rn, self.a[:, None], [1.], node_sizes=(2, 2), num_filters=3),
        ):
            with self.assertRaises(ValueError):
                bad_call()

    def test_dwacd_sro_gating_control_anchors_and_failures(self):
        ctrl = dwacd_sro_gating_control()
        self.assertEqual(ctrl, dwacd_sro_gating_control())
        analytic = ctrl['analytic']
        self.assertAlmostEqual(analytic['true_drift_samples'], 5.12, places=12)
        self.assertAlmostEqual(analytic['true_phase_slope_rad_per_bin'], 0.02 * math.pi, places=14)
        self.assertAlmostEqual(analytic['phase_at_k_min_rad'], 0.32 * math.pi, places=14)
        self.assertAlmostEqual(analytic['phase_at_k_max_rad'], 2.56 * math.pi, places=14)
        self.assertTrue(analytic['wraps_past_pi_at_k_max'])
        self.assertAlmostEqual(analytic['unambiguous_wacd_drift_limit_samples'], 256.0, places=12)
        self.assertAlmostEqual(analytic['intra_frame_drift_samples'], 0.0512, places=14)

        sym = ctrl['symmetric_60pct_speech_40pct_quiet_silence']
        self.assertEqual(sym['jointly_active_pairs'], 60)
        self.assertLess(abs(sym['single_trial']['gated']['error_ppm']), 0.5)
        self.assertLess(sym['gated_rmse_ppm'], 0.5)
        self.assertLess(
            abs(sym['single_trial']['gated']['eps_hat_ppm'] + sym['single_trial']['reversed_conjugation_eps_hat_ppm']),
            1e-5,
        )

        asym = ctrl['asymmetric_60pct_speech_40pct_quiet_silence']
        self.assertEqual(asym['jointly_active_pairs'], 20)
        self.assertAlmostEqual(asym['active_centroid_lag_frames'], 240.0, places=12)
        self.assertAlmostEqual(asym['expected_ungated_centroid_ppm'], 120.0, places=12)
        self.assertLess(abs(asym['single_trial']['gated']['error_ppm']), 0.5)
        self.assertLess(asym['gated_rmse_ppm'], 0.5)
        self.assertAlmostEqual(asym['ungated_mean_bias_ppm'], 20.0, delta=0.5)
        self.assertGreater(asym['ungated_rmse_ppm'], 15.0)

        for key in ('symmetric_60pct_speech_40pct_local_noise_bursts',
                    'sparse_20pct_speech_80pct_local_noise_bursts'):
            case = ctrl[key]
            self.assertLess(abs(case['single_trial']['gated']['error_ppm']), 1.0)
            self.assertLess(case['gated_rmse_ppm'], 0.5)
            self.assertGreater(case['single_trial']['gated']['mean_msc'], 0.8)
            self.assertGreater(
                case['analytic_cross_psd_snr_gated_db'] - case['analytic_cross_psd_snr_ungated_db'],
                25.0,
            )
            self.assertLess(case['single_trial']['ungated']['mean_msc'], 0.05)
            self.assertGreater(case['single_trial']['ungated']['phase_rmse_rad'], 1.5)
            self.assertGreater(case['ungated_rmse_ppm'], 100.0)

        for bad_call in (
            lambda: dwacd_sro_gating_control(seed=True),
            lambda: dwacd_sro_gating_control(trials=0),
            lambda: dwacd_sro_gating_control(trials=True),
        ):
            with self.assertRaises(ValueError):
                bad_call()
        json.dumps(ctrl, allow_nan=False)

    def test_independent_exercise_id_set_and_finite_json(self):
        result = run_experiments()
        expected = {'E13-'+str(i).zfill(2) for i in range(1, 22)}
        self.assertEqual(set(result['exercises']), expected)
        self.assertNotIn('schema_version', result['exercises'])
        json.dumps(result, allow_nan=False)


if __name__ == '__main__':
    unittest.main()
