"""Chapter 1 checks from rational sums, half-angle identities and geometry."""
import json
import math
import hashlib
import io
from pathlib import Path
import subprocess
import sys
import unittest
import tempfile
from unittest.mock import patch
import wave

import numpy as np

from codes.chapters.ch01 import chapter01_experiments as ex
from scripts import run_portable_ci as portable_ci
from codes.chapters.ch00.io_contracts import strict_json_loads, validate_asset_directory
from codes.chapters.ch01.core.spatial_coherence_and_dba import (
    color_white_samples,
    differential_1f_noise_budget,
    equicorrelated_coherence,
    hermitian_psd_sqrt_coloring,
    iec61672_a_weighting,
    sofa_spherical_to_cartesian_and_colatitude,
    spatial_coherence_coloring_and_boundary,
)
from codes.chapters.ch00.core.audio_samples import alignment_error_case, prepare_exports, read_pcm16


class Chapter01ExperimentsTest(unittest.TestCase):
    def test_result_ids_and_finite_json(self):
        if not portable_ci.is_original_generation_host():
            self.skipTest("complete published exercise replay requires the original generation host")
        results = ex.run_exercises()
        self.assertEqual(
            set(results),
            {
                'E01-04', 'E01-05', 'E01-06', 'E01-07', 'E01-08',
                'E01-09', 'E01-10', 'E01-11',
            },
        )
        json.dumps(results, allow_nan=False)

    def test_short_record_cross_terms_are_not_assumed_zero(self):
        full, cropped = ex.finite_record_cross_terms()['cases']
        self.assertEqual(full['products'], [1, -1, -1, 1])
        self.assertEqual(full['average_squared'], [1, 0, 0, 1])
        self.assertEqual(full['output_mean_square'], 1 / 2)
        self.assertEqual(full['cross_second_moment'], 0)
        self.assertEqual(cropped['average_squared'], [1, 0, 0])
        self.assertEqual(cropped['cross_second_moment'], -1 / 3)
        self.assertEqual(cropped['diagonal_contribution'], 1 / 2)
        self.assertEqual(cropped['cross_contribution'], -1 / 6)
        self.assertEqual(cropped['output_mean_square'], 1 / 3)
        self.assertAlmostEqual(cropped['gain_relative_to_channel_1_db'], 10 * math.log10(3))

    def test_residual_delay_has_independent_half_angle_answers(self):
        result = ex.residual_sample_delay()
        low, high = result['cases']
        # cos(pi/16)^2=(2+sqrt(2+sqrt(2)))/4, cos(pi/4)^2=1/2.
        low_power = (2 + math.sqrt(2 + math.sqrt(2))) / 4
        for row, expected in ((low, low_power), (high, .5)):
            self.assertAlmostEqual(row['target_power_ratio'], expected)
            self.assertAlmostEqual(row['projected_target_amplitude_ratio'], math.sqrt(expected), places=12)
            self.assertAlmostEqual(row['target_level_change_db'], 10 * math.log10(expected))
            self.assertEqual(row['aligned_target_amplitude_ratio'], 1)
        self.assertEqual(result['residual_delay_seconds'], 1 / 16000)
        self.assertEqual(result['unaligned_common_delay_samples'], .5)
        self.assertEqual(result['aligned_common_delay_samples'], 1)

    def test_woodworth_zero_side_and_left_right_symmetry(self):
        rows = {row['azimuth_deg']: row for row in ex.woodworth_comparison()['cases']}
        self.assertEqual(rows[0]['sphere_itd_us'], 0)
        self.assertEqual(rows[0]['free_point_itd_us'], 0)
        # At the side, the extra ray length is quarter circumference + radius.
        side_length = .0875 * math.pi / 2 + .0875
        self.assertAlmostEqual(rows[90]['sphere_itd_us'], side_length / 343 * 1e6)
        self.assertAlmostEqual(rows[90]['free_point_itd_us'], .175 / 343 * 1e6)
        self.assertAlmostEqual(rows[30]['free_point_itd_us'], .0875 / 343 * 1e6)
        self.assertAlmostEqual(rows[1]['sphere_itd_us'], 8.90451504010888)
        for angle in (1, 30, 90):
            for key in ('sphere_itd_us', 'free_point_itd_us'):
                self.assertAlmostEqual(rows[-angle][key], -rows[angle][key])

    def test_audio_alignment_is_causal_and_preserves_reference(self):
        case = alignment_error_case()
        signal = case['signals']
        early, later = signal['alignment_array']
        self.assertEqual(later[0], 0)
        np.testing.assert_array_equal(later[1:], early[:-1])
        np.testing.assert_array_equal(signal['alignment_aligned'], later)
        np.testing.assert_array_equal(signal['alignment_reference'], later)
        self.assertEqual(case['parameters']['seed'], None)
        self.assertEqual(signal['alignment_array'].shape, (2, 32000))

    def test_pcm_tone_amplitudes_and_gain_match_independent_identities(self):
        files, groups = prepare_exports({'alignment_error': alignment_error_case()})
        self.assertEqual(len(files), 4)
        self.assertEqual(groups['alignment_error']['common_export_gain'], 1.)
        expected_low = .18 * math.sqrt((2 + math.sqrt(2 + math.sqrt(2))) / 4)
        expected_high = .18 / math.sqrt(2)
        index = np.arange(1600, 30400)
        for stem, expected in [('alignment_reference', [.18, .18]),
                               ('alignment_aligned', [.18, .18]),
                               ('alignment_unaligned', [expected_low, expected_high])]:
            fs, decoded = read_pcm16(files[stem + '.wav'][0])
            self.assertEqual(fs, 16000)
            self.assertEqual(decoded.shape, (1, 32000))
            amplitudes = [2 * abs(np.mean(decoded[0, index] * np.exp(-2j * np.pi * f * index / fs)))
                          for f in (1000, 4000)]
            # DFT magnitude error <= twice the pointwise PCM rounding bound.
            np.testing.assert_allclose(amplitudes, expected, atol=1/32768, rtol=0)
        self.assertEqual(files['alignment_reference.wav'][0], files['alignment_aligned.wav'][0])
        _, stereo = read_pcm16(files['alignment_array.wav'][0])
        self.assertEqual(stereo.shape, (2, 32000))

    def test_published_alignment_scores_match_actual_pcm_real_projection(self):
        result = ex.residual_sample_delay()['published_audio']
        root = Path(__file__).resolve().parents[1]
        index = np.arange(1600, 30400)
        measured = {}
        for stem, record in result['assets'].items():
            raw = (root / record['path']).read_bytes()
            self.assertEqual(hashlib.sha256(raw).hexdigest(), record['sha256'])
            with wave.open(io.BytesIO(raw), 'rb') as stream:
                self.assertEqual(stream.getframerate(), 16000)
                self.assertEqual(stream.getnframes(), 32000)
                self.assertEqual(stream.getsampwidth(), 2)
                channels = stream.getnchannels()
                samples = np.frombuffer(stream.readframes(32000), dtype='<i2').reshape(-1, channels).T / 32768
            # Orthogonal integer-period real projections, rather than the
            # implementation's multi-column least-squares solver.
            amplitudes = []
            for frequency in (1000, 4000):
                phase = 2 * np.pi * frequency * index / 16000
                a = 2 * np.mean(samples[0, index] * np.sin(phase))
                b = 2 * np.mean(samples[0, index] * np.cos(phase))
                amplitudes.append(float(np.hypot(a, b)))
            measured[stem] = amplitudes
            np.testing.assert_allclose(record['tone_amplitudes_by_channel'][0], amplitudes, atol=1e-14, rtol=0)
        float_cases = result['float_measurements']['cases']
        pcm_cases = result['pcm_measurements']['cases']
        expected = [math.sqrt((2 + math.sqrt(2 + math.sqrt(2))) / 4), 1 / math.sqrt(2)]
        for i, (float_case, pcm_case) in enumerate(zip(float_cases, pcm_cases)):
            self.assertAlmostEqual(float_case['unaligned_amplitude_ratio'], expected[i], places=12)
            self.assertAlmostEqual(pcm_case['reference_amplitude'], measured['alignment_reference'][i], places=13)
            self.assertAlmostEqual(pcm_case['unaligned_amplitude_ratio'],
                                   measured['alignment_unaligned'][i] / measured['alignment_reference'][i], places=13)
            self.assertEqual(pcm_case['aligned_amplitude_ratio'], 1.)
        # Quantization has an observable effect; do not silently substitute
        # analytic floating-point values for actual published PCM measurements.
        self.assertGreater(abs(pcm_cases[0]['unaligned_amplitude_ratio'] - expected[0]), 1e-6)
        self.assertEqual(result['scoring_samples'], 28800)

    def test_cross_terms_and_same_mixture_have_independent_rational_answers(self):
        result = ex.mixture_power_and_identifiability()
        for row, cross, total in zip(result['illustrative_cross_terms'], (0, 1, -1), (2, 4, 0)):
            self.assertEqual(row['signal_mean_square'], 1)
            self.assertEqual(row['noise_mean_square'], 1)
            self.assertEqual(row['cross_second_moment'], cross)
            self.assertEqual(row['mixture_mean_square'], total)
            self.assertEqual(row['snr_db'], 0)
        first, second = result['cases']
        self.assertEqual(first['mixture'], second['mixture'])
        self.assertEqual(first['mixture'], [2, 0, -2, 0])
        for row, ps, pv, cross, snr in ((first, 1/2, 1/2, 1/2, 0),
                                       (second, 9/8, 1/8, 3/8, 10 * math.log10(9))):
            self.assertEqual(row['signal_mean_square'], ps)
            self.assertEqual(row['noise_mean_square'], pv)
            self.assertEqual(row['cross_second_moment'], cross)
            self.assertEqual(row['mixture_mean_square'], 2)
            self.assertEqual(ps + pv + 2 * cross, 2)
            self.assertAlmostEqual(row['snr_db'], snr)

    def test_white_noise_gain_counts_actual_noise_and_preserves_target(self):
        equal, aggressive = ex.white_noise_gain()['cases']
        for row, weights, expected_noise, expected_gain in ((equal, [.5, .5], .5, 2),
                                                         (aggressive, [2, -1], 5, .2)):
            self.assertEqual(row['weights'], weights)
            self.assertEqual(row['target_amplitude_gain'], 1)
            self.assertEqual(row['output_noise_mean_square'], expected_noise)
            self.assertEqual(row['white_noise_gain_linear'], expected_gain)
            self.assertAlmostEqual(row['white_noise_gain_db'], 10 * math.log10(expected_gain))

    def test_binaural_exercise_reads_all_five_actual_assets(self):
        from codes.chapters.ch01.examples import generate_binaural_cues as demo
        # A local fixture exercises the actual strict reader in this environment;
        # the catalog test above separately covers published full replay.
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / 'binaural'
            demo.generate_assets(directory)
            with patch.object(ex, 'check_assets', side_effect=lambda: demo.check_assets(directory)):
                result = ex.binaural_level_and_delay()
                wav = directory / 'itd_only.wav'
                changed = bytearray(wav.read_bytes()); changed[244] ^= 1
                wav.write_bytes(changed)
                with self.assertRaises(ValueError):
                    ex.binaural_level_and_delay()
        self.assertEqual(set(result['samples']), {'reference', 'itd_only', 'ild_only', 'consistent', 'conflicting'})
        self.assertEqual(result['common_export_gain'], 1.)
        for name, expected_lag, sign in (('reference', 0, 0), ('itd_only', 8, 0),
                                         ('ild_only', 0, 1), ('consistent', 8, 1), ('conflicting', 8, -1)):
            row = result['samples'][name]
            pcm, floats = row['pcm_measurements'], row['float_measurements']
            self.assertEqual(pcm['estimated_itd_samples'], expected_lag)
            self.assertEqual(pcm['truth_itd_samples'], expected_lag)
            self.assertAlmostEqual(floats['ild_right_minus_left_db'], sign * 20 * math.log10(2))
            self.assertAlmostEqual(pcm['ild_right_minus_left_db'], sign * 20 * math.log10(2), delta=.0002)
            self.assertEqual(pcm['correlation_window_samples'], 28800)
            self.assertEqual(pcm['power_window_samples'], 28800)

    def test_spectral_exercise_reads_four_assets_and_keeps_nonzero_conditions(self):
        from codes.chapters.ch01.examples import generate_spectral_cues as demo
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / 'spectral'
            demo.generate_assets(directory)
            check = demo.check_assets
            with patch.object(demo, 'check_assets', side_effect=lambda: check(directory)):
                result = ex.spectral_response_and_source_spectrum()
                wav = directory / 'flat_stereo.wav'
                changed = bytearray(wav.read_bytes()); changed[244] ^= 1
                wav.write_bytes(changed)
                with self.assertRaises(ValueError):
                    ex.spectral_response_and_source_spectrum()
        manifest = result['published_audio']
        self.assertEqual(set(manifest['files']), {'flat_source.wav', 'flat_stereo.wav',
                                                'tilted_source.wav', 'tilted_stereo.wav'})
        self.assertEqual(manifest['common_export_gain'], 1.)
        self.assertEqual(set(manifest['samples']), {'flat', 'tilted'})
        cosine = math.sqrt(2 + math.sqrt(2)) / 2
        for row, sign in zip(result['analytic_spectral'], (1, -1)):
            self.assertAlmostEqual(row['left_power_gain'], 1.25 + sign*cosine)
            self.assertAlmostEqual(row['right_power_gain'], 1.25 - sign*cosine)
        self.assertTrue(all(row['finite_metrics_rejected'] for row in result['zero_support_controls']))

    def test_spatial_coherence_coloring_matches_rational_and_singular_identities(self):
        result = spatial_coherence_coloring_and_boundary()
        two_real = result['two_mic_real']
        # Independent rational derivation for rho = 24/25:
        # lambda_min = 1/25, lambda_max = 49/25; sqrt = 1/5, 7/5 -> [[4/5, 3/5], [3/5, 4/5]].
        np.testing.assert_allclose(two_real['eigenvalues'], [1 / 25, 49 / 25], atol=1e-14, rtol=0)
        np.testing.assert_allclose(two_real['evd_sqrt_matrix'], [[4 / 5, 3 / 5], [3 / 5, 4 / 5]], atol=1e-14, rtol=0)
        # Cholesky: L21 = 24/25, L22 = sqrt(1 - 576/625) = 7/25.
        np.testing.assert_allclose(two_real['cholesky_lower_matrix'], [[1.0, 0.0], [24 / 25, 7 / 25]], atol=1e-14, rtol=0)
        np.testing.assert_allclose(two_real['evd_colored_samples'],
                                   [[7 / 5, 1 / 5, -1 / 5, -7 / 5], [7 / 5, -1 / 5, 1 / 5, -7 / 5]],
                                   atol=1e-14, rtol=0)
        np.testing.assert_allclose(two_real['cholesky_colored_samples'],
                                   [[1.0, 1.0, -1.0, -1.0], [31 / 25, 17 / 25, -17 / 25, -31 / 25]],
                                   atol=1e-14, rtol=0)
        expected_cov = [[1.0, 24 / 25], [24 / 25, 1.0]]
        np.testing.assert_allclose(two_real['evd_sample_covariance'], expected_cov, atol=1e-14, rtol=0)
        np.testing.assert_allclose(two_real['cholesky_sample_covariance'], expected_cov, atol=1e-14, rtol=0)

        two_complex = result['two_mic_quadrature_complex']
        np.testing.assert_allclose(two_complex['sqrt_real'], [[4 / 5, 0.0], [0.0, 4 / 5]], atol=1e-14, rtol=0)
        np.testing.assert_allclose(two_complex['sqrt_imag'], [[0.0, 3 / 5], [-3 / 5, 0.0]], atol=1e-14, rtol=0)
        np.testing.assert_allclose(two_complex['hermitian_diagonal_real'], [1.0, 1.0], atol=1e-14, rtol=0)
        # Unconjugated transpose gives (4/5)^2 + (3j/5)^2 = 16/25 - 9/25 = 7/25 = 0.28.
        np.testing.assert_allclose(two_complex['unconjugated_transpose_diagonal_real'], [7 / 25, 7 / 25], atol=1e-14, rtol=0)

        four_sing = result['four_mic_singular_boundary']
        self.assertEqual(four_sing['numerical_rank'], 3)
        np.testing.assert_allclose(four_sing['eigenvalues'], [0.0, 4 / 3, 4 / 3, 4 / 3], atol=1e-14, rtol=0)
        diag_val = math.sqrt(3) / 2
        off_val = -math.sqrt(3) / 6
        expected_sqrt_4 = [[diag_val if r == c else off_val for c in range(4)] for r in range(4)]
        np.testing.assert_allclose(four_sing['evd_sqrt_matrix'], expected_sqrt_4, atol=1e-14, rtol=0)
        np.testing.assert_allclose(four_sing['cholesky_pivots_1_to_3'],
                                   [1.0, 2 * math.sqrt(2) / 3, math.sqrt(6) / 3], atol=1e-14, rtol=0)
        self.assertAlmostEqual(four_sing['cholesky_fourth_radicand'], 0.0, places=14)
        np.testing.assert_allclose(four_sing['colored_basis_column_sums'], [0.0, 0.0, 0.0, 0.0], atol=1e-14, rtol=0)
        self.assertTrue(result['non_psd_rejection_control']['rejected'])
        self.assertAlmostEqual(result['non_psd_rejection_control']['min_eigenvalue_analytic'], -1 / 5, places=14)

    def test_spatial_coloring_validates_hermitian_and_psd_contracts(self):
        # Tiny numerical round-off negative eigenvalue within psd_rtol is clamped to 0, while genuine negative raises.
        near_singular = np.array([[1.0, 1.0 + 1e-13], [1.0 + 1e-13, 1.0]], dtype=np.float64)
        clamped = hermitian_psd_sqrt_coloring(near_singular, psd_rtol=1e-10)
        self.assertEqual(clamped['numerical_rank'], 1)
        self.assertEqual(float(clamped['eigenvalues'][0]), 0.0)

        with self.assertRaises(ValueError):
            hermitian_psd_sqrt_coloring(np.array([[1.0, 0.5], [0.2, 1.0]]))
        with self.assertRaises(np.linalg.LinAlgError):
            hermitian_psd_sqrt_coloring(equicorrelated_coherence(4, -0.4))
        with self.assertRaises(ValueError):
            equicorrelated_coherence(True, 0.0)
        with self.assertRaises(ValueError):
            color_white_samples(np.eye(2), np.ones((3, 4)))

    def test_sensor_dba_budget_and_sofa_coordinates_match_closed_forms(self):
        result = ex.sensor_dba_noise_budget()
        spot_by_f = {row['frequency_hz']: row for row in result['iec61672_spot_table']}
        self.assertEqual(spot_by_f[1000.0]['a_db'], 0.0)
        self.assertEqual(spot_by_f[1000.0]['power_weight'], 1.0)
        # Independent direct evaluation of IEC 61672-1 at 100 Hz and 1000 Hz.
        f1, f2, f3, f4 = 20.598997, 107.65265, 737.86223, 12194.217
        ra_100 = (f4**2 * 100.0**4) / ((100.0**2 + f1**2) * math.sqrt((100.0**2 + f2**2) * (100.0**2 + f3**2)) * (100.0**2 + f4**2))
        ra_1k = (f4**2 * 1000.0**4) / ((1000.0**2 + f1**2) * math.sqrt((1000.0**2 + f2**2) * (1000.0**2 + f3**2)) * (1000.0**2 + f4**2))
        self.assertAlmostEqual(result['r_a_1000'], ra_1k, places=14)
        self.assertAlmostEqual(spot_by_f[100.0]['a_db'], 20.0 * math.log10(ra_100 / ra_1k), places=12)
        self.assertAlmostEqual(spot_by_f[100.0]['power_weight'], (ra_100 / ra_1k) ** 2, places=14)

        self.assertEqual(result['mic_b_over_a_unweighted_100hz'], 11 / 2)
        self.assertAlmostEqual(result['mic_b_over_a_unweighted_100hz_db'], 10.0 * math.log10(5.5), places=13)
        self.assertEqual(result['dba_band_50_8000_hz']['num_points'], 796)
        self.assertLess(result['dba_band_50_8000_hz']['a_weighted_ratio_b_over_a'], 1.0)
        self.assertGreater(result['dba_band_50_8000_hz']['apparent_dba_advantage_b_db'], 1.5)
        self.assertEqual(result['dba_band_20_20000_hz']['num_points'], 1999)
        self.assertAlmostEqual(result['dba_band_20_20000_hz']['apparent_dba_advantage_b_db'], 1.9149, places=3)

        diff = result['endfire_differential_1cm']
        expected_wng_100 = 2.0 * math.sin(2.0 * math.pi * 100.0 * 0.01 / 343.0) ** 2
        self.assertAlmostEqual(diff['wng_100hz'], expected_wng_100, places=15)
        # Closed-form asymptotic integral ratio over [100, 4000] Hz is exact rational 49 / 16 = 3.0625.
        self.assertAlmostEqual(diff['speech_band_100_4000_asymptotic_ratio_b_over_a'], 49 / 16, places=14)
        self.assertAlmostEqual(diff['speech_band_100_4000_asymptotic_penalty_b_db'], 20.0 * math.log10(7 / 4), places=13)
        self.assertGreater(diff['speech_band_100_4000_grid_ratio_b_over_a'], 3.1)

        sofa_pt = sofa_spherical_to_cartesian_and_colatitude(-30.0, 30.0, 1.0)
        self.assertEqual(sofa_pt['colatitude_deg'], 60.0)
        self.assertEqual(sofa_pt['book_azimuth_deg'], 30.0)
        self.assertEqual(sofa_pt['azimuth_sofa_wrapped_360_deg'], 330.0)
        np.testing.assert_allclose(sofa_pt['cartesian_sofa_m'], [3 / 4, -math.sqrt(3) / 4, 1 / 2], atol=1e-14, rtol=0)
        np.testing.assert_allclose(sofa_pt['cartesian_book_m'], [math.sqrt(3) / 4, 3 / 4, 1 / 2], atol=1e-14, rtol=0)
        for az_180 in (-180.0, 180.0, 540.0):
            self.assertEqual(
                sofa_spherical_to_cartesian_and_colatitude(az_180, 0.0, 1.0)['book_azimuth_deg'],
                180.0,
            )

        # Descending frequency grid invariance.
        asc_grid = np.array([100.0, 500.0, 1000.0, 4000.0], dtype=np.float64)
        asc_budget = differential_1f_noise_budget(asc_grid, 0.01)
        desc_budget = differential_1f_noise_budget(asc_grid[::-1], 0.01)
        for key in (
            'unweighted_ratio_b_over_a',
            'a_weighted_ratio_b_over_a',
            'differential_grid_ratio_b_over_a',
            'differential_asymptotic_integral_ratio_b_over_a',
            'differential_asymptotic_penalty_b_db',
        ):
            self.assertAlmostEqual(desc_budget[key], asc_budget[key], places=14)
        # Validation error checks on frequency zeros, spatial aliasing, and out-of-range SOFA angles.
        with self.assertRaises(ValueError):
            iec61672_a_weighting(0.0)
        with self.assertRaises(ValueError):
            differential_1f_noise_budget(np.array([17150.0]), 0.01, sound_speed=343.0)
        with self.assertRaises(ValueError):
            sofa_spherical_to_cartesian_and_colatitude(0.0, 95.0, 1.0)
        with self.assertRaises(ValueError):
            sofa_spherical_to_cartesian_and_colatitude(0.0, 0.0, 0.0)


    def _read_published_pcm(self, relative, specification):
        root = Path(__file__).resolve().parents[1]
        directory = root / relative
        validate_asset_directory(directory, set(specification) | {'MANIFEST.json'}, check=True)
        manifest = strict_json_loads((directory / 'MANIFEST.json').read_bytes())
        self.assertEqual(set(manifest['files']), set(specification))
        self.assertEqual(manifest['common_export_gain'], 1.)
        for source, digest in manifest['source_sha256'].items():
            self.assertEqual(hashlib.sha256((root / source).read_bytes()).hexdigest(), digest)
        pcm = {}
        for name, (channels, frames) in specification.items():
            blob = (directory / name).read_bytes()
            self.assertEqual(hashlib.sha256(blob).hexdigest(), manifest['files'][name]['sha256'])
            with wave.open(io.BytesIO(blob), 'rb') as reader:
                self.assertEqual((reader.getframerate(), reader.getnchannels(), reader.getnframes(),
                                  reader.getsampwidth(), reader.getcomptype()),
                                 (16000, channels, frames, 2, 'NONE'))
                raw = reader.readframes(frames)
            self.assertEqual(len(raw), channels * frames * 2)
            pcm[name] = np.frombuffer(raw, dtype='<i2').reshape(-1, channels).T.astype(np.int64)
        return pcm

    def test_published_binaural_pcm_integer_delay_and_gain(self):
        cases = {'reference': ((0, 0), (1, 1)), 'itd_only': ((8, 0), (1, 1)),
                 'ild_only': ((0, 0), (2, 1)), 'consistent': ((8, 0), (2, 1)),
                 'conflicting': ((8, 0), (1, 2))}
        pcm = self._read_published_pcm('codes/chapters/ch01/binaural_audio',
                                      {name + '.wav': (2, 32008) for name in cases})
        reference = pcm['reference.wav'][0, :32000]
        for name, (delays, gain_divisors) in cases.items():
            data = pcm[name + '.wav']
            energies = []
            for channel, (delay, divisor) in enumerate(zip(delays, gain_divisors)):
                self.assertTrue(np.all(data[channel, :delay] == 0))
                self.assertTrue(np.all(data[channel, delay + 32000:] == 0))
                # Two independently rounded copies of x and x/2 differ by at
                # most one integer count after doubling the latter.
                error = divisor * data[channel, delay:delay + 32000] - reference
                self.assertLessEqual(int(np.max(abs(error))), 0 if divisor == 1 else 1)
                segment = data[channel, 1600 + delay:30400 + delay]
                energies.append(int(segment @ segment))
            expected_ild = 20 * math.log10(gain_divisors[0] / gain_divisors[1])
            self.assertAlmostEqual(10 * math.log10(energies[1] / energies[0]), expected_ild, delta=.0002)
            right = data[1, 1600:30400]
            scores = []
            for lag in range(-16, 17):
                left = data[0, 1600 + lag:30400 + lag]
                cross, left_energy, right_energy = int(left @ right), int(left @ left), int(right @ right)
                scores.append(cross / math.sqrt(left_energy * right_energy))
            self.assertEqual(int(np.argmax(scores)) - 16, delays[0] - delays[1])

    def test_published_spectral_pcm_integer_fir_bound(self):
        pcm = self._read_published_pcm('codes/chapters/ch01/spectral_audio',
                                      {name + suffix + '.wav': (channels, frames)
                                       for name in ('flat', 'tilted')
                                       for suffix, channels, frames in (('_source', 1, 32000), ('_stereo', 2, 32001))})
        for name in ('flat', 'tilted'):
            source = pcm[name + '_source.wav'][0]
            current = np.pad(source, (0, 1)); previous = np.pad(source, (1, 0))
            stereo = pcm[name + '_stereo.wav']
            for channel, sign in ((0, 1), (1, -1)):
                # |round(y) - round(x) - .5*round(x_prev)| <= 1.25 counts.
                # Doubling gives an integer residual with magnitude <= 2.
                residual = 2 * stereo[channel] - 2 * current - sign * previous
                self.assertLessEqual(int(np.max(abs(residual))), 2)
                self.assertGreater(int(stereo[channel] @ stereo[channel]), 0)

    def test_direct_script_bootstrap_with_published_assets(self):
        if not portable_ci.is_original_generation_host():
            self.skipTest("direct script evaluates all published assets on the original generation host")
        script = Path(__file__).resolve().parents[1] / 'codes/chapters/ch01/chapter01_experiments.py'
        proc = subprocess.run([sys.executable, str(script)], capture_output=True, text=True, check=True)
        self.assertEqual(set(json.loads(proc.stdout)), set(ex.run_exercises()))


if __name__ == '__main__':
    unittest.main()
