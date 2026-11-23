"""Analytic/fractional expectations and actual PCM for E08-10..27."""
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import wave
import numpy as np
from codes.chapters.ch08.chapter08_experiments import (
    ROOT,
    auxiva_iss_vs_ip_experiment,
    css_published_audio,
    run_experiments,
)
from codes.chapters.ch08.core.separation import si_sdr
from codes.chapters.ch08.examples import generate_css_audio, mask_representation_demo
from scripts.run_portable_ci import is_original_generation_host
from codes.chapters.ch08.core.separation_metrics_streaming import (
    ci_sdr_toeplitz,
    single_speaker_bypass_gate,
    streaming_subband_mac_budget,
)


def complex_array(value):
    return np.array(value['real']) + 1j * np.array(value['imag'])


class Chapter08Exercises(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if is_original_generation_host():
            cls.r = run_experiments()
            return
        with tempfile.TemporaryDirectory() as temporary:
            mask = Path(temporary) / 'mask'
            css = Path(temporary) / 'css'
            mask_representation_demo.generate_assets(mask)
            generate_css_audio.generate_assets(css)
            check_mask = mask_representation_demo.check_assets
            check_css = generate_css_audio.check_assets
            with patch.object(mask_representation_demo, 'check_assets', side_effect=lambda: check_mask(mask)), \
                 patch.object(generate_css_audio, 'check_assets', side_effect=lambda: check_css(css)):
                cls.r = run_experiments()

    def test_ids_and_finite_json(self):
        self.assertEqual(set(self.r), {f'E08-{i:02}' for i in range(10, 28)})
        json.dumps(self.r, allow_nan=False)

    def test_sequential_complex_rows(self):
        r = self.r['E08-10']
        expected = [[np.sqrt(2/3), -1j/np.sqrt(6)], [0, 1/np.sqrt(2)]]
        np.testing.assert_allclose(complex_array(r['sequential_rows']), expected, atol=1e-15)
        np.testing.assert_allclose(complex_array(r['simultaneous_weighted_gram']), [[1, -.5j], [.5j, 1]], atol=1e-15)

    def test_joint_vs_scalar_regression(self):
        r = self.r['E08-11']
        np.testing.assert_allclose(r['joint_gains'], [1, 2])
        np.testing.assert_allclose(r['individual_gains'], [2, 5/2])
        np.testing.assert_allclose(r['individual_sum'], [2, 5/2, 9/2])
        np.testing.assert_allclose(complex_array(r['complex_application_gain']), [2+1j])

    def test_new_constraint_order_and_background_controls(self):
        r=self.r['E08-22']
        self.assertEqual(r['order_max_waveform_difference'],1/8)
        self.assertEqual(r['consistency_then_mix']['consistency_max_spectrum_error'],1/4)
        self.assertEqual(r['mix_then_consistency']['consistency_max_spectrum_error'],0)
        b=self.r['E08-23']['background_constraint_control']
        np.testing.assert_allclose(complex_array(b['J']),[[4/5]])
        np.testing.assert_allclose(complex_array(b['background_second_moment']),[[42/25]])

    def test_css_continuous_post_overlap_gain_control(self):
        r=self.r['E08-24']
        self.assertEqual(r['overlap_matching']['current_indices_for_previous'],[0,1])
        np.testing.assert_allclose(r['overlap_application_gain'],[-.5,-.5])
        for name,nmse in (('reference',0),('naive',9),('polarity',1),('corrected',0)):
            a=r['analytic'][name]['post_overlap']
            np.testing.assert_allclose(a['reference_NMSE_per_channel'],[nmse,nmse],atol=1e-12)
        self.assertEqual(len(r['published_audio']['files']),4)
        # Actual PCM corrected/reference identity is stronger than rounded
        # display scores and is checked without regenerating either waveform.
        directory=ROOT/'codes/chapters/ch08/css_audio'
        self.assertEqual((directory/'css_corrected.wav').read_bytes(),
                         (directory/'css_reference.wav').read_bytes())

    def test_conditional_prior_counterexample(self):
        r = self.r['E08-12']
        np.testing.assert_allclose(r['posterior'], [[9/10, 1/10], [0, 1], [0, 1]])
        np.testing.assert_allclose(r['mean_posterior_priors'], [.3, .7])
        self.assertAlmostEqual(r['change'], np.log(17/25))
        self.assertLess(r['change'], 0)

    def test_nmf_substep_objectives(self):
        r = self.r['E08-13']
        np.testing.assert_allclose(r['updated_basis'], [[1], [3]])
        np.testing.assert_allclose(r['updated_activation'], [[np.sqrt(2), np.sqrt(2)]])
        np.testing.assert_allclose(r['objectives_initial_basis_activation'], [20, 8+2*np.log(3), 4*np.sqrt(2)+2*np.log(2)+2*np.log(3)])

    def test_underdetermined_covariance(self):
        r = self.r['E08-14']
        np.testing.assert_allclose(r['observation_scm_a'], np.eye(2))
        np.testing.assert_allclose(r['observation_scm_b'], np.eye(2))
        np.testing.assert_allclose(r['first_wiener_image_a'], [.5, 1])
        np.testing.assert_allclose(r['first_wiener_image_b'], [.6, .4])

    def test_congruence_wiener_fractions(self):
        r = self.r['E08-15']
        np.testing.assert_allclose(r['commutator'], [[0, -5], [5, 0]])
        np.testing.assert_allclose(r['microphone_images'], [[11/12, 2/3], [13/12, 1/3]])
        np.testing.assert_allclose(r['sum_images'], [2, 1])

    def test_bounded_real_masks(self):
        r = self.r['E08-21']['unit_scale_complex_control']
        np.testing.assert_allclose(complex_array(r['complex_oracle_mask']), [5, (1+1j)/2])
        np.testing.assert_allclose(r['squared_errors'], [16/25, 1/2])

    def test_overlap_matching_and_ambiguities(self):
        r = self.r['E08-16']
        m = r['audio_parameters']['matching']
        np.testing.assert_allclose(m['absolute_centered_correlation'], [[20/101, 1], [1, 20/101]], atol=1e-13)
        self.assertEqual(m['current_indices_for_previous'], [1, 0])
        for key in ('silence', 'identical_slots'):
            self.assertEqual(r[key]['status'], 'ambiguous')
            self.assertIsNone(r[key]['current_indices_for_previous'])

    def test_published_css_pcm_score_and_sources(self):
        r = self.r['E08-16']['published_audio']
        self.assertEqual(r['score_interval_samples'], [0, 32000])
        self.assertEqual(r['common_export_gain'], 1)
        measured = r['actual_pcm_measurements']
        expected = {'css_overlap_naive': [-3.340758788514897, -3.3405346277825743],
                    'css_overlap_aligned': [19.998887508407464, 19.999652603596694]}
        for stem, values in expected.items():
            scores = measured[stem]['fixed_identity_scores']
            np.testing.assert_allclose([item['si_sdr_db'] for item in scores], values, atol=1e-11)
            self.assertEqual(measured[stem]['samples_per_channel'], 32000)
            for score in scores:
                self.assertIsInstance(score['reference_centered_energy_numerator'], int)
                self.assertGreater(score['residual_energy_ratio_denominator'], 0)
        for name, digest in r['source_sha256'].items():
            self.assertEqual(hashlib.sha256((ROOT/name).read_bytes()).hexdigest(), digest)

    def test_css_stale_source_and_corrupt_asset_are_not_repaired(self):
        original = ROOT/'codes/chapters/ch00/audio/MANIFEST.json'
        audio = ROOT/'codes/chapters/ch08/audio'
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            manifest = json.loads(original.read_text())
            name = 'codes/chapters/ch08/core/css.py'
            manifest['generator_inputs'][name] = '0'*64
            path = directory/'MANIFEST.json'
            path.write_text(json.dumps(manifest))
            before = path.read_bytes()
            with self.assertRaisesRegex(ValueError, 'source SHA'):
                css_published_audio(manifest_path=path, audio_directory=audio)
            self.assertEqual(path.read_bytes(), before)
            path.write_bytes(original.read_bytes())
            for stem in ('reference', 'mixture', 'naive', 'aligned'):
                filename = f'css_overlap_{stem}.wav'
                (directory/filename).write_bytes((audio/filename).read_bytes())
            corrupted = directory/'css_overlap_aligned.wav'
            blob = bytearray(corrupted.read_bytes())
            blob[-1] ^= 1
            corrupted.write_bytes(blob)
            with self.assertRaisesRegex(ValueError, 'metadata or SHA'):
                css_published_audio(manifest_path=path, audio_directory=directory)
            self.assertEqual(corrupted.read_bytes(), bytes(blob))

    def test_css_nonfinite_boolean_or_incomplete_scores_are_rejected(self):
        original = ROOT/'codes/chapters/ch00/audio/MANIFEST.json'
        audio = ROOT/'codes/chapters/ch08/audio'
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'MANIFEST.json'
            for scores in ([float('nan'), 0.], [float('inf'), 0.],
                           [True, 0.], [0.], ['0', 0.]):
                with self.subTest(scores=scores):
                    manifest = json.loads(original.read_text())
                    manifest['groups']['css_overlap']['parameters']['pcm_si_sdr_db'][
                        'css_overlap_naive'] = scores
                    path.write_text(json.dumps(manifest))
                    before = path.read_bytes()
                    with self.assertRaises(ValueError):
                        css_published_audio(manifest_path=path, audio_directory=audio)
                    self.assertEqual(path.read_bytes(), before)

    def test_whitened_four_points_are_dependent(self):
        r = self.r['E08-17']
        np.testing.assert_allclose(r['mixture_covariance'], [[5/4, 1], [1, 5/4]])
        np.testing.assert_allclose(r['whitened_points'], r['source_points'], atol=1e-15)
        np.testing.assert_allclose(r['rotated_covariance'], np.eye(2), atol=1e-15)
        self.assertEqual(r['probability_first_zero'], 1/2)
        self.assertEqual(r['probability_second_zero'], 1/2)
        self.assertEqual(r['probability_both_zero'], 0)
        self.assertEqual(r['independent_product_probability'], 1/4)
        self.assertEqual(np.linalg.matrix_rank(r['rank_one_covariance']), 1)

    def test_three_spd_joint_congruence_obstruction(self):
        r = self.r['E08-18']
        np.testing.assert_array_equal(r['whitened_remaining_commutator'], [[0, -1], [1, 0]])
        np.testing.assert_allclose(r['shape_eigenvalues'], [[1, 1], [1, 2], [1, 3]])
        self.assertFalse(r['joint_congruence_exists_for_all_three'])

    def test_fastmnmf_jacobian_cancels_model_scale(self):
        first, second = self.r['E08-19']['cases']
        expected = 7/12+np.log(12)
        for case in (first, second):
            self.assertAlmostEqual(case['quadratic_term'], 7/12)
            self.assertAlmostEqual(case['full_objective'], expected)
            np.testing.assert_allclose(case['physical_covariance'], [[7, 3], [3, 3]])
        self.assertAlmostEqual(second['jacobian_term'], -2*np.log(4))
        self.assertAlmostEqual(second['objective_without_jacobian']-first['objective_without_jacobian'], 2*np.log(4))

    def test_weighted_projection_kkt_fractions(self):
        r = self.r['E08-20']
        np.testing.assert_array_equal(r['weighted_projection'], [3/2, 5/2])
        np.testing.assert_array_equal(r['equal_projection'], [2, 2])
        np.testing.assert_array_equal(r['kkt_correction_over_variance'], [1/2, 1/2])
        self.assertEqual(r['weighted_cost_at_weighted_projection'], 1/2)
        self.assertAlmostEqual(r['weighted_cost_at_equal_projection'], 2/3)

    def test_mask_representation_reads_published_pcm(self):
        directory = ROOT / 'codes/chapters/ch08/mask_audio'
        manifest = json.loads((directory / 'MANIFEST.json').read_text(encoding='utf-8'))
        self.assertEqual(manifest['sample_rate_hz'], 16000)
        self.assertEqual(manifest['common_export_gain'], 1)
        self.assertEqual(len(manifest['files']), 6)
        samples = {}
        for filename, record in manifest['files'].items():
            path = directory / filename
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), record['sha256'])
            with wave.open(str(path), 'rb') as stream:
                self.assertEqual((stream.getnchannels(), stream.getframerate(),
                                  stream.getsampwidth(), stream.getnframes()), (1, 16000, 2, 32000))
                pcm = np.frombuffer(stream.readframes(32000), dtype='<i2').astype(np.int64)
            samples[filename] = pcm[2400:29600]

        records = manifest['samples']
        reference = samples[records['target']['file']]
        self.assertEqual(int(reference @ reference), 292062085900)
        for stem, numerator in [('bounded_real', 166479798100),
                                ('unbounded_real', 73012418550),
                                ('complex_oracle', 0)]:
            sample = records[stem]
            error = samples[sample['file']] - reference
            self.assertEqual(int(error @ error), numerator)
            measured = sample['pcm_measurements']
            self.assertEqual(measured['integer_error_squared_sum'], numerator)
            self.assertEqual(measured['integer_reference_squared_sum'], 292062085900)
            self.assertEqual(measured['integer_mse_denominator'], 27200*32768**2)

    def test_matched_slots_can_cancel_and_have_gain_error(self):
        r = self.r['E08-24']['four_point_discrete_control']
        p = np.array([[1, -1, 0, 0], [0, 0, 1, -1]])
        first, second = r['cases']
        for case in (first, second):
            self.assertEqual(case['matching']['status'], 'matched')
            self.assertEqual(case['matching']['current_indices_for_previous'], [0, 1])
            np.testing.assert_array_equal(case['signed_centered_correlation'], -np.eye(2))
            np.testing.assert_array_equal(case['corrected_equal_overlap'], p)
        np.testing.assert_array_equal(first['naive_equal_overlap'], np.zeros((2, 4)))
        np.testing.assert_array_equal(first['direct_application_gain'], [-1, -1])
        np.testing.assert_array_equal(second['naive_equal_overlap'], -.5*p)
        np.testing.assert_array_equal(second['polarity_only_overlap'], 1.5*p)
        np.testing.assert_array_equal(second['direct_application_gain'], [-.5, -.5])

    def test_ci_sdr_vs_si_sdr_and_overfitting_boundary(self):
        r = self.r['E08-25']['cases']
        # Case 1: pure 1-sample causal delay s_1 = [0, 1, 0, -1] against s_0 = [1, 0, -1, 0].
        delay_case = r['pure_one_sample_delay']
        self.assertEqual(delay_case['si_sdr_L1']['ci_sdr_db'], -120.0)
        self.assertEqual(delay_case['ci_sdr_L2']['ci_sdr_db'], 120.0)
        np.testing.assert_allclose(delay_case['ci_sdr_L2']['filter_taps'], [0.0, 1.0], atol=1e-15)
        np.testing.assert_array_equal(delay_case['ci_sdr_L2']['gram_matrix'], [[2.0, 0.0], [0.0, 2.0]])

        # Case 2: 2-tap room response s_0 + s_1 + w, where ||s_0||^2=2, ||s_1||^2=2, ||w||^2=1/25.
        room_case = r['two_tap_room_plus_distortion']
        expected_l1_ratio = F(2, 1) / (F(2, 1) + F(1, 25))  # 50/51
        self.assertEqual(expected_l1_ratio, F(50, 51))
        self.assertAlmostEqual(
            room_case['si_sdr_L1']['ci_sdr_db'],
            10.0 * math.log10(float(expected_l1_ratio)),
            places=12,
        )
        np.testing.assert_allclose(room_case['ci_sdr_L2']['filter_taps'], [1.0, 1.0], atol=1e-15)
        self.assertAlmostEqual(room_case['ci_sdr_L2']['target_energy'], 4.0, places=14)
        self.assertAlmostEqual(room_case['ci_sdr_L2']['residual_energy'], 0.04, places=14)
        self.assertAlmostEqual(room_case['ci_sdr_L2']['ci_sdr_db'], 20.0, places=12)

        # Case 3: a cropped square fit overfits, but the full convolution tail
        # adds independent residual energy. Exact normal equations give thirds.
        overfit_case = r['full_tail_vs_truncated_L_equals_N']
        self.assertEqual(overfit_case['truncated_L4_determinant'], 1.0)
        self.assertAlmostEqual(overfit_case['gram_L4_determinant'], 9.0)
        np.testing.assert_allclose(overfit_case['truncated_L4_filter_taps'], [1, -1, 2, -2])
        self.assertEqual(overfit_case['truncated_L4_residual_energy'], 0.0)
        self.assertEqual(overfit_case['si_sdr_L1']['ci_sdr_db'], -120.0)
        self.assertEqual(overfit_case['ci_sdr_L2']['ci_sdr_db'], -120.0)
        self.assertFalse(overfit_case['ci_sdr_L4']['overparameterized'])
        self.assertEqual(overfit_case['ci_sdr_L4']['rank'], 4)
        np.testing.assert_allclose(overfit_case['ci_sdr_L4']['filter_taps'],
                                   [float(F(x, 3)) for x in (1, -1, 2, -2)], atol=1e-14)
        self.assertAlmostEqual(overfit_case['ci_sdr_L4']['target_energy'], float(F(4, 3)))
        self.assertAlmostEqual(overfit_case['ci_sdr_L4']['residual_energy'], float(F(8, 3)))
        self.assertAlmostEqual(overfit_case['ci_sdr_L4']['ci_sdr_db'], 10 * math.log10(.5))

        # Verify L=1 agreement with Chapter 8 si_sdr and strict invalid-input rejection.
        vec_est = np.array([1.2, -0.4, -0.8, 0.5])
        vec_ref = np.array([1.0, -0.5, -1.0, 0.5])
        self.assertAlmostEqual(
            ci_sdr_toeplitz(vec_est, vec_ref, filter_length=1, zero_mean=True)['ci_sdr_db'],
            si_sdr(vec_est, vec_ref, zero_mean=True),
            places=12,
        )
        for bad_args, bad_kwargs in (
            ((np.zeros(4), vec_ref), {'filter_length': 1}),
            ((vec_est, np.ones(4)), {'filter_length': 1, 'zero_mean': True}),
            ((vec_est, vec_ref), {'filter_length': 0}),
            ((vec_est, vec_ref), {'filter_length': 5}),
            ((vec_est, vec_ref), {'filter_length': True}),
        ):
            with self.assertRaises(ValueError):
                ci_sdr_toeplitz(*bad_args, **bad_kwargs)

    def test_streaming_mac_budget_and_bypass_gate(self):
        r = self.r['E08-26']
        mac = r['mac_budget']
        # Independent exact integer MAC/s: (16000 // 128) * 4 * F * 4 * 64^2.
        self.assertEqual(mac['fullband_256']['macs_per_second'], 125 * 4 * 256 * 4 * 4096)
        self.assertEqual(mac['fullband_256']['macs_per_second'], 2_097_152_000)
        self.assertEqual(mac['mel_band_64']['macs_per_second'], 125 * 4 * 64 * 4 * 4096)
        self.assertEqual(mac['mel_band_64']['macs_per_second'], 524_288_000)
        self.assertEqual(mac['compression_factor'], 4.0)

        fric = r['fricative_splitting_case']
        self.assertAlmostEqual(fric['slot1_uncentered_si_sdr_db'], 10.0 * math.log10(99.0), places=12)
        self.assertAlmostEqual(
            fric['bypass_gate']['secondary_to_primary_db'],
            10.0 * math.log10(1.0 / 99.0),
            places=12,
        )
        self.assertTrue(fric['bypass_gate']['bypass'])
        self.assertAlmostEqual(fric['bypassed_reconstruction_error_energy'], 0.0, places=14)
        self.assertEqual(fric['unbypassed_wer'], float(F(1, 4)))
        self.assertEqual(fric['bypassed_wer'], 0.0)

        ovl = r['two_speaker_overlap_case']['bypass_gate']
        self.assertAlmostEqual(ovl['secondary_to_primary_db'], -10.0, places=12)
        self.assertFalse(ovl['bypass'])

        # Exact -18 dB boundary, primary_idx == 1 preservation, and invalid inputs.
        at_boundary = np.array([[1.0, 0.0], [0.0, 10.0 ** (-0.9)]])
        self.assertTrue(single_speaker_bypass_gate(at_boundary, threshold_db=-18.0)['bypass'])
        above_boundary = np.array([[1.0, 0.0], [0.0, 10.0 ** (-0.89)]])
        self.assertFalse(single_speaker_bypass_gate(above_boundary, threshold_db=-18.0)['bypass'])
        slot1_primary = np.array([[0.0, 1.0, 0.0, 0.0], [math.sqrt(99.0), 0.0, 0.0, 0.0]])
        gate_slot1 = single_speaker_bypass_gate(slot1_primary, threshold_db=-18.0)
        self.assertEqual(gate_slot1['primary_slot_index'], 1)
        self.assertEqual(gate_slot1['secondary_slot_index'], 0)
        self.assertTrue(gate_slot1['bypass'])
        np.testing.assert_array_equal(gate_slot1['output_waveforms'][0], np.zeros(4))
        np.testing.assert_allclose(
            gate_slot1['output_waveforms'][1],
            slot1_primary[0] + slot1_primary[1],
            atol=1e-15,
        )
        for bad_call in (
            lambda: single_speaker_bypass_gate(np.zeros((2, 4))),
            lambda: single_speaker_bypass_gate(at_boundary, threshold_db=0.0),
            lambda: streaming_subband_mac_budget(16000, 0, 64, 4, 64),
        ):
            with self.assertRaises(ValueError):
                bad_call()

    def test_auxiva_iss_rank1_vs_ip_convergence(self):
        r = self.r['E08-27']
        self.assertEqual(r['sweeps'], 5)
        self.assertEqual(r['initial_cost'], 8.0)
        self.assertEqual(r['ip_matrix_solves_per_sweep'], 2)
        self.assertEqual(r['iss_matrix_solves_per_sweep'], 0)

        # Sweep 1 IP exact analytic checks:
        ip1, ip2 = r['sweep1_ip']
        np.testing.assert_allclose(
            complex_array(ip1['unnormalized_direction']),
            [float(F(5, 16)), 1j * float(F(1, 16))],
            atol=1e-15,
        )
        self.assertAlmostEqual(ip1['cost_after'], 5.0 + math.log(float(F(16, 5))), places=13)
        np.testing.assert_allclose(
            complex_array(ip2['unnormalized_direction']),
            [-1j * float(F(1, 20)), float(F(3, 10))],
            atol=1e-15,
        )
        self.assertAlmostEqual(ip2['cost_after'], 2.0 + math.log(float(F(320, 29))), places=13)
        expected_w_ip1 = np.array([
            [math.sqrt(5.0) / 4.0, -1j / (4.0 * math.sqrt(5.0))],
            [1j / (2.0 * math.sqrt(29.0)), 3.0 / math.sqrt(29.0)],
        ])
        np.testing.assert_allclose(complex_array(ip2['W_after']), expected_w_ip1, atol=1e-14)

        # Sweep 1 ISS exact analytic checks (zero matrix inversions):
        iss1, iss2 = r['sweep1_iss']
        np.testing.assert_allclose(
            complex_array(iss1['steering_vector']),
            [float(F(1, 2)), -1j * float(F(1, 5))],
            atol=1e-15,
        )
        np.testing.assert_allclose(
            complex_array(iss1['W_after']),
            [[0.5, 0.0], [0.2j, 1.0]],
            atol=1e-15,
        )
        self.assertAlmostEqual(iss1['cost_after'], float(F(21, 5)) + 2.0 * math.log(2.0), places=13)
        np.testing.assert_allclose(
            complex_array(iss2['steering_vector']),
            [1j * float(F(5, 58)), 1.0 - math.sqrt(5.0) / 4.0],
            atol=1e-14,
        )
        expected_w_iss1 = np.array([
            [float(F(15, 29)), -1j * float(F(5, 58))],
            [1j / (4.0 * math.sqrt(5.0)), math.sqrt(5.0) / 4.0],
        ])
        np.testing.assert_allclose(complex_array(iss2['W_after']), expected_w_iss1, atol=1e-14)
        self.assertAlmostEqual(
            iss2['cost_after'],
            float(F(54, 29)) + math.log(float(F(64, 5))),
            places=13,
        )

        # Monotone descent and 5-sweep convergence to exact global minimum Q*.
        self.assertTrue(r['ip_monotone_decreasing'])
        self.assertTrue(r['iss_monotone_decreasing'])
        self.assertGreater(r['initial_cost'], r['iss_cost_history'][1])
        self.assertGreater(r['initial_cost'], r['ip_cost_history'][1])
        q_star = 2.0 + math.log(64.0 * (3.0 - 2.0 * math.sqrt(2.0)))
        self.assertAlmostEqual(r['exact_global_min_cost'], q_star, places=14)
        self.assertLess(r['final_cost_gap_ip_vs_iss'], 1e-10)
        self.assertLess(r['ip_gap_to_optimum'], 1e-10)
        self.assertLess(r['iss_gap_to_optimum'], 1e-10)

        # Final weighted Grams diagonalize (V1, V2) with unit target diagonals.
        mu2 = 3.0 + 2.0 * math.sqrt(2.0)
        for key in ('ip_final_weighted_grams', 'iss_final_weighted_grams'):
            g1, g2 = (complex_array(item) for item in r[key])
            np.testing.assert_allclose(g1, np.diag([1.0, mu2]), atol=1e-6)
            np.testing.assert_allclose(g2, np.diag([mu2, 1.0]), atol=1e-6)

        for bad_sweeps in (0, -1, True, 2.5):
            with self.assertRaises(ValueError):
                auxiva_iss_vs_ip_experiment(sweeps=bad_sweeps)  # type: ignore[arg-type]


if __name__ == '__main__':
    unittest.main()
