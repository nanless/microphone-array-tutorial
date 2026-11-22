"""Independent state, circular support and retarded-time regression anchors."""
import copy
import tempfile
from pathlib import Path
import unittest
import numpy as np

from codes.chapters.ch09.core.tracking import (
    ConstantVelocityKalman, CircularParticleFilter, _checked_covariance, systematic_resample, wrap_angle,
)
from codes.chapters.ch09.core.moving_source import retarded_emission_times, synthetic_source, free_field_array
from codes.chapters.ch00.core.audio_samples import pcm16_bytes
from codes.chapters.ch00.cross_chapter.enhancement_structure_exercises import imm_mix
from codes.chapters.ch09.examples.moving_source_audio import build_fixture, generate


class TrackingBoundaries(unittest.TestCase):
    def test_negative_variance_is_never_tolerated(self):
        for scale in (1e-200, 1., 1e200):
            for field in ('covariance', 'process_noise'):
                arguments = {'covariance': np.eye(2)*scale, 'process_noise': np.zeros((2, 2))}
                arguments[field] = np.diag([-1e-13, 1.])*scale
                with self.assertRaisesRegex(ValueError, 'negative variance'):
                    ConstantVelocityKalman([0., 0.], **arguments)

    def test_relative_symmetry_and_significant_negative_direction(self):
        for scale in (1e-200, 1., 1e200):
            for matrix in (np.array([[1., .001], [0., 1.]]), np.array([[1., 1.01], [1.01, 1.]])):
                with self.assertRaises(ValueError):
                    _checked_covariance(matrix*scale)

    def test_singular_psd_and_smallest_subnormal_preserved(self):
        for scale in (1e-200, 1., 1e200):
            repaired = _checked_covariance(np.ones((3, 3))*scale)
            np.testing.assert_allclose(repaired/scale, np.ones((3, 3)), atol=2e-14)
            self.assertGreaterEqual(np.linalg.eigvalsh(repaired/scale).min(), -1e-14)
        tiny = np.nextafter(0., 1.)
        np.testing.assert_array_equal(_checked_covariance(np.diag([tiny, tiny])), np.diag([tiny, tiny]))
        np.testing.assert_array_equal(_checked_covariance(np.zeros((2, 2))), np.zeros((2, 2)))

    def test_joseph_retains_small_observation_variance(self):
        # Scalar posterior is P*R/(P+R), which is 1 to float precision.
        tracker = ConstantVelocityKalman([0., 0.], np.diag([1e20, 1.]), np.zeros((2, 2)))
        tracker.update(0., 1.)
        self.assertAlmostEqual(tracker.covariance[0, 0], 1.)
        self.assertEqual(1e20-(1e20/(1e20+1))*1e20, 0.)

    def test_large_initial_angle_is_same_direction(self):
        tracker = ConstantVelocityKalman([1e20, 0.], np.eye(2), np.zeros((2, 2)))
        self.assertEqual(tracker.state[0], -80.)
        tracker.update(0., 1.)
        np.testing.assert_allclose(tracker.state, [-40., 0.])

    def test_kalman_failure_does_not_commit(self):
        tracker = ConstantVelocityKalman([1., 1e308], np.eye(2), np.zeros((2, 2)))
        before = tracker.state.copy(), tracker.covariance.copy()
        with self.assertRaises(ValueError):
            tracker.predict(2.)
        np.testing.assert_array_equal(tracker.state, before[0])
        np.testing.assert_array_equal(tracker.covariance, before[1])

    def test_particle_large_observation_uses_canonical_direction(self):
        for angle in (1e20, 1e100, 1e308):
            first, second = CircularParticleFilter([0., 90.]), CircularParticleFilter([0., 90.])
            first.update(angle, 10., clutter_probability=0.)
            second.update(wrap_angle(angle), 10., clutter_probability=0.)
            np.testing.assert_array_equal(first.weights, second.weights)
        # For -80 deg, errors are 80 and 170: odds exp(-(170²-80²)/200).
        np.testing.assert_allclose(CircularParticleFilter([179., -179.]).estimate(), -180.)

    def test_positive_prior_support_survives_narrow_likelihood(self):
        pf = CircularParticleFilter([0., 10.])
        pf.update(10., .01, clutter_probability=0.)
        np.testing.assert_array_equal(pf.weights, [0., 1.])
        self.assertEqual(pf.weight_underflow_indices, [0])
        self.assertTrue(np.isfinite(pf.log_weights[0]))
        pf.update(0., .01, clutter_probability=0.)
        np.testing.assert_allclose(pf.weights, [.5, .5], atol=1e-10)
        for std in (1.5e-154, 1e-300, np.nextafter(0., 1.)):
            before = pf.weights.copy(), pf.log_weights.copy()
            with self.assertRaisesRegex(ValueError, 'log likelihood'):
                pf.update(0., std, clutter_probability=0.)
            np.testing.assert_array_equal(pf.weights, before[0])
            np.testing.assert_array_equal(pf.log_weights, before[1])

    def test_particle_small_and_large_std_are_finite(self):
        for std in (np.nextafter(0., 1.), 1e-300, 1e308):
            for clutter in (0., .05, np.nextafter(0., 1.)):
                pf = CircularParticleFilter([0., 10.])
                if std < 1e-150 and clutter == 0:
                    with self.assertRaises(ValueError):
                        pf.update(0., std, clutter_probability=clutter)
                    continue
                with np.errstate(all='raise'):
                    pf.update(0., std, clutter_probability=clutter)
                self.assertTrue(np.all(np.isfinite(pf.weights)))
                self.assertAlmostEqual(pf.weights.sum(), 1.)

    def test_particle_prediction_failure_preserves_rng_and_state(self):
        for velocity, dt, std in ((1e308, 2., 0.), (0., 1., 1e308)):
            pf = CircularParticleFilter(np.arange(100.))
            rng = np.random.default_rng(0)
            state = copy.deepcopy(rng.bit_generator.state)
            before = pf.particles.copy()
            with self.assertRaises(ValueError):
                pf.predict(velocity, dt, std, rng)
            self.assertEqual(rng.bit_generator.state, state)
            np.testing.assert_array_equal(pf.particles, before)

    def test_successful_prediction_consumes_same_random_draws(self):
        pf = CircularParticleFilter([0., 90.]); rng = np.random.default_rng(4); independent = np.random.default_rng(4)
        expected = wrap_angle(np.array([0., 90.])+2.+independent.normal(0., .5, 2))
        pf.predict(20., .1, .5, rng)
        np.testing.assert_array_equal(pf.particles, expected)
        self.assertEqual(rng.bit_generator.state, independent.bit_generator.state)

    def test_systematic_resample_upper_boundary_and_trailing_zero_weights(self):
        class FixedRNG:
            def __init__(self, u: float):
                self._u = u
            def random(self) -> float:
                return self._u

        u_max = float(np.nextafter(1.0, 0.0))
        # (u_max + n - 1) / n rounds to 1.0 in float64 for n >= 2 (e.g., n = 3, 7, 10).
        for n in (1, 2, 3, 7, 10, 100):
            weights = np.ones(n, dtype=float)
            indices = systematic_resample(weights, FixedRNG(u_max))
            self.assertEqual(indices.shape, (n,))
            self.assertTrue(np.all((indices >= 0) & (indices < n)))
            if n >= 2:
                pf = CircularParticleFilter(np.linspace(-90.0, 90.0, n))
                prior = np.zeros(n, dtype=float)
                prior[0], prior[1] = 0.8, 0.2
                pf.set_prior_weights(prior)
                self.assertTrue(pf.resample_if_needed(FixedRNG(u_max), threshold=n))
                self.assertEqual(pf.particles.shape, (n,))
                self.assertTrue(set(pf.last_resample_ancestor_indices.tolist()).issubset({0, 1}))

        # Trailing zero weights must never be chosen, even when cumsum(positive) < 1.0 in float64.
        trailing_cases = (
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([1.0] * 9 + [0.0, 0.0]),  # sum(1/9 * 9) = 0.9999999999999999 < 1.0
            np.array([0.1, 0.2, 0.3, 0.4, 0.0, 0.0]),
            np.array([0.0, 0.25, 0.0, 0.75, 0.0]),
        )
        for weights in trailing_cases:
            positive_indices = set(np.flatnonzero(weights > 0.0).tolist())
            for u in (0.0, 0.5, u_max):
                indices = systematic_resample(weights, FixedRNG(u))
                self.assertEqual(indices.shape, weights.shape)
                self.assertTrue(set(indices.tolist()).issubset(positive_indices))

    def test_unnormalized_mutated_weights_normalize_log_weights_and_ess(self):
        pf = CircularParticleFilter([0.0, 30.0, 60.0])
        pf.weights = np.array([2.0, 6.0, 0.0])
        expected_probs = np.array([0.25, 0.75, 0.0])
        np.testing.assert_allclose(pf.log_weights[:2], np.log(expected_probs[:2]), atol=1e-15)
        self.assertTrue(np.isneginf(pf.log_weights[2]))
        self.assertAlmostEqual(float(np.sum(np.exp(pf.log_weights))), 1.0, places=15)
        self.assertAlmostEqual(pf.effective_sample_size, 1.0 / (0.25**2 + 0.75**2), places=15)

        # In-place mutation scaling weights by a constant leaves normalized log_weights and ESS invariant.
        pf.weights[:] = [10.0, 30.0, 0.0]
        np.testing.assert_allclose(pf.log_weights[:2], np.log(expected_probs[:2]), atol=1e-15)
        self.assertAlmostEqual(pf.effective_sample_size, 1.6, places=15)

        # Huge finite unnormalized weights avoid overflow in sum.
        pf.weights = np.array([1e308, 1e308, 0.0])
        np.testing.assert_allclose(pf.log_weights[:2], [-np.log(2.0), -np.log(2.0)], atol=1e-15)
        self.assertAlmostEqual(pf.effective_sample_size, 2.0, places=15)

        # Zero-sum and negative mutated weights are rejected by both log_weights and effective_sample_size.
        for bad in (np.array([0.0, 0.0, 0.0]), np.array([-1.0, 2.0, 0.0]), np.array([np.nan, 1.0, 0.0])):
            pf.weights = bad
            with self.assertRaises(ValueError):
                _ = pf.log_weights
            with self.assertRaises(ValueError):
                _ = pf.effective_sample_size

    def test_imm_legal_singular_covariance_and_real_contract(self):
        for scale in (1e-200, 1., 1e200):
            out = imm_mix([1.], [[1.]], [[0., 0., 0.]], [np.ones((3, 3))*scale])
            np.testing.assert_allclose(out['covariances'][0]/scale, np.ones((3, 3)), atol=2e-14)
        for invalid in ([True], ['1'], [1+0j]):
            with self.assertRaises(ValueError):
                imm_mix(invalid, [[1]], [[0]], [[[1]]])


class MovingBoundaryTests(unittest.TestCase):
    def test_analytic_quadratic_motion_root(self):
        times = np.array([0., .2, 1., 2.]); p0 = np.array([-.8, 1.5]); velocity = np.array([.8, 0.]); mics = np.array([[-.05, 0.], [.05, 0.]])
        # Solve flight time s=t-u from a quadratic, independently of Newton.
        q = p0 + times[None, :, None]*velocity - mics[:, None, :]
        dot = q@velocity; a = 343.**2-velocity@velocity
        travel = (-dot+np.sqrt(dot**2+a*np.sum(q*q, axis=2)))/a
        expected = times-travel
        actual = retarded_emission_times(times, mics, source_start_xy=p0, source_velocity_xy=velocity)
        np.testing.assert_allclose(actual, expected, atol=2e-15, rtol=0)

    def test_nonfinite_intermediates_rejected_not_nan_success(self):
        with self.assertRaises(ValueError):
            retarded_emission_times([0., 1.], [[0., 0.], [.1, 0.]], source_start_xy=[1e200, 1.], source_velocity_xy=[0., 0.])
        with self.assertRaises(ValueError):
            synthetic_source([1e308])

    def test_nonreal_inputs_and_scalar_arrays_rejected(self):
        for bad in (np.array([0.+2j]), np.array([True]), np.array(['0'])):
            with self.assertRaises(ValueError):
                retarded_emission_times(bad, [[0., 0.], [.1, 0.]], source_start_xy=[0., 1.], source_velocity_xy=[0., 0.])
        for speed in (True, [343.], 343+0j):
            with self.assertRaises(ValueError):
                retarded_emission_times([0.], [[0., 0.], [.1, 0.]], source_start_xy=[0., 1.], source_velocity_xy=[0., 0.], sound_speed=speed)
        with self.assertRaises(ValueError):
            free_field_array([0.], [[0., 0.], [.1, 0.]], source_start_xy=[0., 1.], source_velocity_xy=[0., 0.], reference_distance_m=True)

    def test_existing_three_wav_bytes_preserved(self):
        signals, metadata = build_fixture()
        directory = Path(__file__).resolve().parents[1]/'codes/chapters/ch09/moving_audio'
        for stem, waveform in signals.items():
            self.assertEqual(pcm16_bytes(waveform, 16000), (directory/f'{stem}.wav').read_bytes())

    def test_moving_source_read_only_check_detects_staleness(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            generate(directory)
            generate(directory, check=True)
            path = directory/'source.wav'
            path.write_bytes(path.read_bytes()[:-1])
            before = {p.name: p.read_bytes() for p in directory.iterdir()}
            with self.assertRaises(ValueError):
                generate(directory, check=True)
            self.assertEqual(before, {p.name: p.read_bytes() for p in directory.iterdir()})


if __name__ == '__main__':
    unittest.main()
