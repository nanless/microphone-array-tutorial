"""Chapter 9 anchors (E09-10..22) independently derived from fractions/moments/geometry."""
from fractions import Fraction as F
import math
from pathlib import Path
import struct
import tempfile
import unittest
from unittest.mock import patch
import wave
import numpy as np
from codes.chapters.ch09.chapter09_experiments import (
    actual_tracking_audio, run_experiments, small_set_distances,
)
from codes.chapters.ch09.core.echo_av_fusion import (
    _wrap_angle_deg,
    clear_mot_angular_summary,
    echo_sector_gate,
    fov_bernoulli_missed_detection_update,
    oosm_cv_retrodiction_update,
)
from codes.chapters.ch09.examples.chapter09_tracking_audio import generate


class Chapter09ExperimentsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.audio_directory = Path(cls.temporary.name)/'tracking_audio'
        # Only this temporary fixture is generated; published assets are read-only.
        generate(cls.audio_directory)
        cls.results = run_experiments(audio_directory=cls.audio_directory)

    def test_thirteen_stable_ids(self):
        self.assertEqual(set(self.results), {f'E09-{i:02}' for i in range(10, 23)})

    def test_bearing_ekf_geometry_and_fraction_posterior(self):
        item = self.results['E09-10']
        self.assertAlmostEqual(item['bearing_rad'], math.atan(.5))
        np.testing.assert_allclose(item['jacobian_rad_per_m'], [.4, -.2])
        np.testing.assert_allclose(item['gain_m_per_rad'], [5/3, -5/6])
        np.testing.assert_allclose(item['updated_position_m'], [7/6, 23/12])
        np.testing.assert_allclose(item['posterior_covariance_m2'], [[1/12, 1/12], [1/12, 5/24]])
        # Independent central differences for the atan2(x,y) observation.
        position = np.array([1., 2.]); eps = 1e-5
        derivatives = []
        for direction in np.eye(2):
            high, low = position+eps*direction, position-eps*direction
            derivatives.append((math.atan2(*high)-math.atan2(*low))/(2*eps))
        np.testing.assert_allclose(derivatives, [.4, -.2], atol=1e-10)

    def test_ut_matches_gaussian_second_and_fourth_moments(self):
        item = self.results['E09-11']
        # For a standard normal E[X^2]=1, E[X^4]=3; odd moments vanish.
        self.assertEqual(item['mean'], 1.)
        self.assertEqual(item['variance'], 3-1**2)
        self.assertEqual(item['state_observation_cross_covariance'], 0.)
        self.assertEqual(item['covariance_weights'], [2., .5, .5])

    def test_circular_mean_and_supported_particle_boundary(self):
        item = self.results['E09-12']
        self.assertEqual(item['circular_mean_deg'], -180.)
        self.assertAlmostEqual(item['resultant_length'], math.cos(math.pi/180))
        self.assertIsNone(item['antipodal_angle_deg'])
        self.assertEqual(item['antipodal_status'], 'undefined_antipodal_mean')
        self.assertEqual(item['support_after_first'], [0., 1.])
        self.assertEqual(item['linear_underflow_indices_after_first'], [0])
        np.testing.assert_allclose(item['log_weights_after_first'], [-500000., 0.], atol=1e-8)
        np.testing.assert_allclose(item['support_after_second'], [.5, .5], atol=1e-10)
        np.testing.assert_allclose(item['log_weights_after_second'], [-math.log(2)]*2, atol=1e-10)
        self.assertEqual(item['explicit_zero_prior_after_update'], [0., 1.])
        self.assertEqual(item['prior_reset_method'], 'set_prior_weights')

    def test_phd_cardinality_and_bernoulli_empty_observation(self):
        item = self.results['E09-13']
        for key, variance in [('model_A_count_probability_0_1_2', 0), ('model_B_count_probability_0_1_2', 1)]:
            probabilities = item[key]
            mean = sum(i*p for i, p in enumerate(probabilities))
            self.assertEqual(mean, 1)
            self.assertEqual(sum((i-mean)**2*p for i, p in enumerate(probabilities)), variance)
        np.testing.assert_allclose(item['empty_observation_posterior_mass'], [.1, .1])
        bernoulli = item['bernoulli_empty_observation']
        r = F(4, 5)
        for case, detection, posterior in zip(bernoulli['cases'], [F(9, 10), F(1, 10)], [F(2, 7), F(18, 23)]):
            numerator = r*(1-detection)
            denominator = (1-r)+numerator
            self.assertEqual(numerator/denominator, posterior)
            self.assertAlmostEqual(case['empty_observation_probability'], float(denominator))
            self.assertAlmostEqual(case['bernoulli_posterior_existence'], float(posterior))
            self.assertAlmostEqual(case['different_poisson_phd_posterior_mass'], float(numerator))

    def test_ospa_and_gospa_cardinality_normalization_and_boundaries(self):
        item = self.results['E09-14']
        self.assertEqual(item['original']['ospa'], 6.)
        self.assertEqual(item['original']['gospa_alpha2'], 7.)
        self.assertEqual(item['added_common_correct_target']['ospa'], 4.)
        self.assertEqual(item['added_common_correct_target']['gospa_alpha2'], 7.)
        self.assertEqual(small_set_distances([], [30, 70])['ospa'], 10.)
        self.assertEqual(small_set_distances([179], [-179])['ospa'], 2.)
        self.assertAlmostEqual(small_set_distances([32], [30, 70], order=2)['ospa'], math.sqrt(52))
        self.assertAlmostEqual(small_set_distances([32], [30, 70], order=2)['gospa_alpha2'], math.sqrt(54))
        self.assertEqual(small_set_distances([30, 70], [32]), small_set_distances([32], [30, 70]))
        for first, second, kwargs in [([0, 360], [0], {}), ([0], [1], {'order': .5}), ([0], [1], {'cutoff': 0})]:
            with self.assertRaises(ValueError):
                small_set_distances(first, second, **kwargs)

    def test_fixed_lag_rts_and_oosm_late_observation_joint_gaussian(self):
        item = self.results['E09-15']
        # RTS part: x1 prior variance 2, x2 variance 3, cov 2; measurement covariance [[3,2],[2,4]].
        cross = np.array([2., 2.]); innovation_covariance = np.array([[3., 2.], [2., 4.]])
        posterior_mean = cross@np.linalg.solve(innovation_covariance, [1., 0.])
        posterior_variance = 2-cross@np.linalg.solve(innovation_covariance, cross)
        self.assertAlmostEqual(posterior_mean, .5)
        self.assertAlmostEqual(posterior_variance, .5)
        self.assertAlmostEqual(item['smoothed_time1_mean'], posterior_mean)
        self.assertAlmostEqual(item['smoothed_time1_variance'], posterior_variance)
        np.testing.assert_allclose(item['filtered_means'], [2/3, 1/4])
        np.testing.assert_allclose(item['filtered_variances'], [2/3, 5/8])
        self.assertEqual(item['extra_observation_delay_steps'], 1)
        # OOSM part: independent inverse of S=[[3,2],[2,4]], determinant 8.
        oosm = item['oosm_late_observation']
        inverse = [[F(1, 2), F(-1, 4)], [F(-1, 4), F(3, 8)]]
        oosm_cross = [F(2), F(3)]
        gain = [sum(oosm_cross[i]*inverse[i][j] for i in range(2)) for j in range(2)]
        self.assertEqual(gain, [F(1, 4), F(5, 8)])
        mean = 2*gain[0]-gain[1]
        variance = 3-sum(g*x for g, x in zip(gain, oosm_cross))
        self.assertEqual((mean, variance), (F(-1, 8), F(5, 8)))
        self.assertAlmostEqual(oosm['joint_mean'], float(mean))
        self.assertAlmostEqual(oosm['joint_variance'], float(variance))
        self.assertAlmostEqual(oosm['chronological_replay'][-1]['mean'], float(mean))
        self.assertAlmostEqual(oosm['chronological_replay'][-1]['variance'], float(variance))
        self.assertEqual(oosm['drop_late_observation'], {'mean': -.75, 'variance': .75})
        self.assertAlmostEqual(oosm['wrong_current_time_update']['mean'], 3/7)
        self.assertAlmostEqual(oosm['wrong_current_time_update']['variance'], 3/7)

    def test_audio_uses_independent_denominators_and_measured_timing(self):
        item = self.results['E09-16']
        self.assertEqual(item['pcm_scores']['valid_observation_count'], 173)
        self.assertEqual(item['pcm_scores']['missing_observation_count'], 24)
        self.assertAlmostEqual(item['first_state_time_s'], 255.5/16000)
        self.assertEqual(item['first_available_time_s'], 512/16000)
        self.assertAlmostEqual(item['first_available_time_s']-item['first_state_time_s'], 256.5/16000)
        self.assertNotEqual(item['pcm_scores']['raw_valid_rmse_deg'], item['float_scores']['raw_valid_rmse_deg'])
        self.assertEqual(item['asset_check'], 'strict_read_only_passed')

    def test_audio_scores_decode_actual_file_not_replay_pcm(self):
        path = self.audio_directory/'array_noisy.wav'
        with wave.open(str(path), 'rb') as reader:
            count, channels = reader.getnframes(), reader.getnchannels()
            data = reader.readframes(count)
        integers = struct.unpack('<'+str(count*channels)+'h', data)
        expected = np.array([[integers[i*channels+j]/32768 for i in range(count)]
                             for j in range(channels)])
        from codes.chapters.ch09.core.tracking_audio import analyze_array
        with patch('codes.chapters.ch09.chapter09_experiments.analyze_array', wraps=analyze_array) as score:
            actual_tracking_audio(self.audio_directory)
        score.assert_called_once()
        np.testing.assert_array_equal(score.call_args.args[0], expected)

    def test_audio_tampering_is_rejected_without_repair(self):
        for member in ('array_noisy.wav', 'MANIFEST.json'):
            with self.subTest(member=member), tempfile.TemporaryDirectory() as root:
                directory = Path(root)/'tracking_audio'
                generate(directory)
                path = directory/member
                changed = bytearray(path.read_bytes()); changed[-1] ^= 1
                path.write_bytes(changed)
                before = {p.name: p.read_bytes() for p in directory.iterdir()}
                with self.assertRaises(ValueError):
                    actual_tracking_audio(directory)
                self.assertEqual(before, {p.name: p.read_bytes() for p in directory.iterdir()})

    def test_small_set_tiny_positive_cutoff_and_power_status(self):
        cutoff = 1e-200
        empty = small_set_distances([], [30], cutoff=cutoff, order=2)
        self.assertEqual(empty['ospa'], cutoff)
        self.assertAlmostEqual(empty['gospa_alpha2']/cutoff, 1/math.sqrt(2))
        near = small_set_distances([0], [cutoff/2], cutoff=cutoff, order=2)
        self.assertEqual(near['ospa'], cutoff/2)
        self.assertEqual(near['gospa_alpha2'], cutoff/2)
        self.assertIsNone(near['matched_capped_cost_power'])
        self.assertEqual(near['matched_capped_cost_power_status'], 'positive_power_underflow')
        exact_zero = small_set_distances([0], [0], cutoff=cutoff, order=2)
        self.assertEqual(exact_zero['matched_capped_cost_power'], 0)
        self.assertEqual(exact_zero['matched_capped_cost_power_status'], 'representable')
        with self.assertRaises(ValueError):
            small_set_distances([0], [math.ulp(0.)], cutoff=180, order=2)

    def test_correlated_observations_and_bearing_scale_nullspace(self):
        item = self.results['E09-17']
        corr = item['correlated_observations']
        rho = F(9, 10)
        determinant = 4-(1+rho)**2
        gain = (1-rho)/determinant
        self.assertEqual(gain, F(10, 39))
        np.testing.assert_allclose(corr['joint_gain'], [float(gain)]*2)
        self.assertAlmostEqual(corr['joint_mean'], float(2*gain))
        self.assertAlmostEqual(corr['joint_variance'], float(1-2*gain))
        self.assertAlmostEqual(corr['wrong_independent_mean'], 2/3)
        self.assertAlmostEqual(corr['wrong_independent_variance'], 1/3)
        self.assertGreater(corr['joint_variance'], corr['wrong_independent_variance'])

        scale = item['bearing_scale_ambiguity']
        np.testing.assert_array_equal(scale['bearings_deg'], scale['scaled_bearings_deg'])
        np.testing.assert_allclose(scale['ranges_m'], [math.sqrt(5), math.sqrt(8), math.sqrt(13)])
        np.testing.assert_array_equal(scale['scaled_ranges_m'], 2*np.array(scale['ranges_m']))
        rows = [[2, -1, 0, 0], [2, -2, 2, -2], [2, -3, 4, -6]]
        np.testing.assert_array_equal(scale['stacked_jacobian_without_positive_row_denominators'], rows)
        a, b, c = [[F(row[j]) for j in (0, 1, 3)] for row in rows]
        det3 = a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
        self.assertEqual(det3, 4)
        self.assertEqual(scale['stacked_rank'], 3)
        self.assertTrue(all(sum(F(v)*u for v, u in zip(row, [1, 2, 1, 0])) == 0 for row in rows))
        np.testing.assert_array_equal(scale['stacked_null_product'], [0, 0, 0])

    def test_imm_lifecycle_and_spherical_vmf_present(self):
        imm = self.results['E09-18']
        self.assertEqual(len(imm['steps']), 3)
        self.assertAlmostEqual(imm['steps'][0]['fused_variance'], 0.897358874054, places=9)
        self.assertAlmostEqual(imm['steps'][1]['fused_variance'], 2.556447667640, places=9)
        self.assertAlmostEqual(imm['steps'][2]['fused_variance'], 4.117809823150, places=9)
        life = self.results['E09-19']
        self.assertEqual(life['lifecycle']['counts']['confirmed_publications'], 187)
        self.assertEqual(life['state_age_comparison']['counts']['confirmed_publications'], 189)
        sph = self.results['E09-20']
        self.assertAlmostEqual(sph['world_angles']['azimuth_deg'], -60.0, places=12)
        self.assertAlmostEqual(sph['static_vmf_posterior']['concentration'], 5.0, places=12)

    def test_echo_sector_gating_hijack_and_double_talk(self):
        item = self.results['E09-21']
        hijack = item['ungated_hijack']
        # Step 1 exact fractions: P1^- = [[6, 1], [1, 1]], S1 = 10, K1 = [3/5, 1/10], nu1 = -65.
        self.assertEqual(hijack['step1_nis'], float(F((-65) ** 2, 10)))
        self.assertTrue(hijack['nis_gate_9_rejects_step1'])
        np.testing.assert_allclose(hijack['step1_posterior_state'], [-19.0, -6.5], atol=1e-14)
        np.testing.assert_allclose(
            hijack['step1_posterior_covariance'],
            [[float(F(12, 5)), float(F(2, 5))], [float(F(2, 5)), float(F(9, 10))]],
            atol=1e-14,
        )
        # Step 2 exact fractions: P2^- = [[51/10, 13/10], [13/10, 9/10]], S2 = 91/10, K2 = [51/91, 13/91], nu2 = -39/2.
        s2_angle = F(-51, 2) + F(51, 91) * F(-39, 2)
        s2_vel = F(-13, 2) + F(13, 91) * F(-39, 2)
        self.assertEqual((s2_angle, s2_vel), (F(-255, 7), F(-65, 7)))
        np.testing.assert_allclose(
            hijack['step2_posterior_state'],
            [float(F(-255, 7)), float(F(-65, 7))],
            atol=1e-13,
        )
        self.assertEqual(hijack['unmasked_3hit_false_track_births'], 1)

        gated = item['sector_gated_single_talk']
        self.assertTrue(gated['gate_active'])
        self.assertEqual(gated['retained_detections_deg'], [])
        self.assertEqual(gated['suppressed_detections_deg'], [-45.0])
        np.testing.assert_allclose(gated['step1_coasted_covariance'], [[6.0, 1.0], [1.0, 1.0]])
        np.testing.assert_allclose(gated['step2_coasted_covariance'], [[10.0, 2.0], [2.0, 1.0]])
        self.assertEqual(gated['false_track_births'], 0)

        dt_ctrl = item['double_talk_control']
        self.assertFalse(dt_ctrl['gate_active'])
        np.testing.assert_allclose(dt_ctrl['near_end_22deg_posterior_state'], [float(F(106, 5)), float(F(1, 5))])

        # Wrap-around sector boundary check near -180/+180 deg.
        wrapped_gate = echo_sector_gate(
            np.array([176.0, 160.0]),
            speaker_angle_deg=-178.0,
            sector_half_width_deg=8.0,
            far_end_active=True,
            double_talk=False,
        )
        np.testing.assert_allclose(wrapped_gate['retained_detections_deg'], [160.0])
        np.testing.assert_allclose(wrapped_gate['suppressed_detections_deg'], [176.0])
        with self.assertRaises(ValueError):
            echo_sector_gate(
                np.array([0.0]),
                speaker_angle_deg=0.0,
                sector_half_width_deg=10.0,
                far_end_active=1,  # type: ignore[arg-type]
                double_talk=False,
            )

    def test_av_oosm_retrodiction_fov_and_mota_vs_ospa(self):
        item = self.results['E09-22']
        oosm = item['oosm_retrodiction']
        naive = oosm['naive_current_time']
        self.assertAlmostEqual(naive['innovation_deg'], -3.0, places=12)
        np.testing.assert_allclose(naive['posterior_state'], [float(F(103, 10)), 94.0], atol=1e-13)

        retro = oosm['retrodiction']
        self.assertAlmostEqual(retro['innovation_deg'], 0.0, places=14)
        self.assertAlmostEqual(retro['innovation_variance_deg2'], float(F(889, 100)), places=12)
        np.testing.assert_allclose(
            retro['kalman_gain'],
            [float(F(840, 889)), float(F(1700, 889))],
            atol=1e-13,
        )
        np.testing.assert_allclose(retro['posterior_state'], [13.0, 100.0], atol=1e-13)
        p11 = F(9) - F(840, 889) * F(42, 5)
        p12 = F(20) - F(840, 889) * F(17)
        p22 = F(100) - F(1700, 889) * F(17)
        self.assertEqual((p11, p12, p22), (F(135, 127), F(500, 127), F(60000, 889)))
        expected_cov = [[float(p11), float(p12)], [float(p12), float(p22)]]
        np.testing.assert_allclose(retro['posterior_covariance'], expected_cov, atol=1e-12)
        np.testing.assert_allclose(oosm['chronological_replay_state_at_t4'], [13.0, 100.0], atol=1e-13)
        np.testing.assert_allclose(oosm['chronological_replay_covariance_at_t4'], expected_cov, atol=1e-12)

        fov = item['camera_fov_existence']
        self.assertFalse(fov['in_fov'])
        self.assertAlmostEqual(fov['fov_aware_posterior_existence'], float(F(4, 5)), places=14)
        naive_r = F(4, 5) * F(1, 20) / (1 - F(4, 5) * F(19, 20))
        self.assertEqual(naive_r, F(1, 6))
        self.assertAlmostEqual(fov['naive_posterior_existence'], float(naive_r), places=14)

        swap = item['identity_vs_set_metrics']['pure_id_swap']
        self.assertEqual(swap['ospa_per_frame_deg'], [0.0, 0.0])
        self.assertEqual(swap['clear_mot']['idsw'], 2)
        self.assertEqual(swap['clear_mot']['fn'], 0)
        self.assertEqual(swap['clear_mot']['fp'], 0)
        self.assertAlmostEqual(swap['clear_mot']['mota'], float(F(1, 2)), places=14)

        swap_echo = item['identity_vs_set_metrics']['id_swap_plus_echo_false_track']
        self.assertAlmostEqual(swap_echo['ospa_per_frame_deg'][1], float(F(10, 3)), places=12)
        self.assertEqual(swap_echo['clear_mot']['fp'], 1)
        self.assertEqual(swap_echo['clear_mot']['idsw'], 2)
        self.assertAlmostEqual(swap_echo['clear_mot']['mota'], float(F(1, 4)), places=14)

        reappear = clear_mot_angular_summary(
            [{1: 0.0}, {2: 0.0}, {1: 2.0, 2: 0.0}],
            [{1: 0.0}, {1: 0.0}, {1: 0.0, 2: 2.0}],
        )
        self.assertEqual(reappear['idsw'], 0)
        self.assertEqual(reappear['motp_deg'], 0.0)
        self.assertEqual(reappear['mota'], 1.0)

        self.assertEqual(_wrap_angle_deg(180.0), -180.0)
        self.assertEqual(_wrap_angle_deg(-180.0), -180.0)
        wrapped_oosm = oosm_cv_retrodiction_update(
            np.array([-170.0, 200.0]), np.eye(2), 170.0, 1.0, 0.1
        )
        self.assertAlmostEqual(
            wrapped_oosm['retrodiction']['predicted_measurement_deg'], 170.0, places=12
        )

        fov_outside_unit = fov_bernoulli_missed_detection_update(
            1.0, 60.0, fov_half_width_deg=45.0, in_fov_detection_prob=1.0
        )
        self.assertFalse(fov_outside_unit['in_fov'])
        self.assertEqual(fov_outside_unit['fov_aware_posterior_existence'], 1.0)
        self.assertIsNone(fov_outside_unit['naive_posterior_existence'])
        fov_offset_cam = fov_bernoulli_missed_detection_update(
            0.8, 100.0, fov_half_width_deg=20.0, in_fov_detection_prob=0.5, camera_center_deg=90.0
        )
        self.assertTrue(fov_offset_cam['in_fov'])
        self.assertAlmostEqual(fov_offset_cam['fov_aware_posterior_existence'], float(F(2, 3)), places=14)

        for bad_call in (
            lambda: oosm_cv_retrodiction_update(np.array([0.0, 0.0]), np.eye(2), 0.0, 1.0, -0.01),
            lambda: oosm_cv_retrodiction_update(
                np.array([0.0, 0.0]), np.array([[1.0, 2.0], [2.0, 1.0]]), 0.0, 1.0, 0.1
            ),
            lambda: fov_bernoulli_missed_detection_update(1.0, 0.0, in_fov_detection_prob=1.0),
            lambda: clear_mot_angular_summary([{}], [{}]),
        ):
            with self.assertRaises(ValueError):
                bad_call()

    def test_oosm_retrodiction_accepts_rank1_psd_with_roundoff(self):
        # Construct rank-1 PSD covariances whose floating-point eigvalsh can dip slightly below zero.
        for vec in (
            np.array([1.0, 3.0]),
            np.array([1.0, 1.0 / 3.0]),
            np.array([3.0, 7.0]),
            np.array([1e3, 3e3]),
        ):
            cov = np.outer(vec, vec)
            res = oosm_cv_retrodiction_update(np.array([10.0, 20.0]), cov, 12.0, 1.5, 0.1)
            for branch in ('retrodiction', 'naive_current_time'):
                post = res[branch]['posterior_covariance']
                np.testing.assert_allclose(post, post.T, atol=1e-15, rtol=0)
                eigvals = np.linalg.eigvalsh(post)
                self.assertGreaterEqual(eigvals[0], -1e-14 * max(1.0, float(eigvals[-1])))

        # Explicit tiny negative roundoff within tolerance is projected cleanly; genuine indefiniteness fails.
        eps = np.finfo(float).eps
        slightly_neg = np.array([[1.0, 1.0], [1.0, 1.0 - 4.0 * eps]])
        res = oosm_cv_retrodiction_update(np.array([0.0, 0.0]), slightly_neg, 1.0, 1.0, 0.1)
        self.assertGreaterEqual(np.linalg.eigvalsh(res['retrodiction']['posterior_covariance'])[0], -1e-14)
        with self.assertRaises(ValueError):
            oosm_cv_retrodiction_update(
                np.array([0.0, 0.0]), np.array([[1.0, 1.0], [1.0, 1.0 - 1e-8]]), 1.0, 1.0, 0.1
            )


if __name__ == '__main__':
    unittest.main()
