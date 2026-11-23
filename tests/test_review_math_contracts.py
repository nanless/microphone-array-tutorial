"""Independent controls for the mathematical contracts found in the book review."""

from fractions import Fraction
import math
import unittest
import warnings

import numpy as np

from codes.chapters.ch00.core import chapter00_math as math_core
from codes.chapters.ch02.core.conventions import hermitian_part
from codes.chapters.ch04.core.covariance import _load_covariance, _validate_covariance


class ReviewedCovarianceScales(unittest.TestCase):
    def test_already_hermitian_subnormal_matrices_keep_their_values(self):
        unit = np.nextafter(0.0, 1.0)
        for count in (1, 2, 3, 5):
            matrix = count * unit * np.array([[2, 1j], [-1j, 2]])
            with self.subTest(count=count):
                np.testing.assert_array_equal(hermitian_part(matrix), matrix)
                np.testing.assert_array_equal(_validate_covariance(matrix), matrix)

    def test_component_midpoints_match_exact_binary_rationals_at_extremes(self):
        unit, largest = np.nextafter(0.0, 1.0), np.finfo(float).max
        values = (0.0, unit, 3 * unit, 5 * unit, np.finfo(float).tiny,
                  1.0, np.nextafter(1.0, 2.0), largest, -largest, -3 * unit)
        for first in values:
            for second in values:
                expected = float((Fraction.from_float(first)
                                  + Fraction.from_float(second)) / 2)
                matrix = np.array([[0, complex(first, first)],
                                   [complex(second, -second), 0]])
                actual = hermitian_part(matrix)
                self.assertEqual(actual[0, 1], complex(expected, expected))
                self.assertEqual(actual[1, 0], complex(expected, -expected))

    def test_relative_loading_retains_small_positive_trace_and_rounding(self):
        unit = np.nextafter(0.0, 1.0)
        for diagonal, mean_units in (([1, 1, 1], 1), ([2, 2, 2], 2),
                                     ([1, 1, 1, 7], 2), ([1, 1, 1, 11], 4)):
            matrix = np.diag(np.asarray(diagonal) * unit)
            loaded = _load_covariance(matrix, 1.0, 1e12)
            expected = np.diag((np.asarray(diagonal) + mean_units) * unit)
            np.testing.assert_array_equal(loaded, expected)

    def test_finite_large_mean_is_available_when_the_trace_overflows(self):
        level = np.finfo(float).max / 2
        matrix = level * np.eye(3)
        with warnings.catch_warnings():
            warnings.simplefilter('error', RuntimeWarning)
            loaded = _load_covariance(matrix, 0.25, 1e12)
        np.testing.assert_allclose(loaded / level, 1.25 * np.eye(3), rtol=3e-16)


class ReviewedAnalyticalComplexMoments(unittest.TestCase):
    def test_propriety_uses_analytic_parameters_across_scales(self):
        for scale in (1.0, 1e-20, 1e200):
            self.assertTrue(math_core.complex_second_order_moments(scale, scale)['is_proper_second_order'])
            self.assertFalse(math_core.complex_second_order_moments(scale, 0)['is_proper_second_order'])
            self.assertFalse(math_core.complex_second_order_moments(scale, 2 * scale)['is_proper_second_order'])
            self.assertFalse(math_core.complex_second_order_moments(scale, scale, 1e-16)['is_proper_second_order'])
        self.assertTrue(math_core.complex_second_order_moments(0, 0, 1)['is_proper_second_order'])

    def test_large_and_mixed_variances_keep_representable_moments(self):
        for rho in (0.0, 0.5):
            result = math_core.complex_second_order_moments(1e200, 1e200, rho)
            self.assertEqual(result['variance_r'], 2e200)
            self.assertTrue(math.isfinite(result['pseudo_variance_abs']))
            self.assertAlmostEqual(result['pseudo_variance_c'].imag / 1e200, 2 * rho)
        result = math_core.complex_second_order_moments(1e-300, 1e300, 0.5)
        self.assertAlmostEqual(result['pseudo_variance_c'].imag, 1.0)

    def test_unrepresentable_nonzero_moments_do_not_return_false_zeros(self):
        unit = np.nextafter(0.0, 1.0)
        with self.assertRaisesRegex(ValueError, 'component'):
            math_core.complex_second_order_moments(unit, unit, 0.1)
        with self.assertRaisesRegex(ValueError, 'total variance'):
            math_core.complex_second_order_moments(np.finfo(float).max, np.finfo(float).max)

    def test_subnormal_geometric_mean_is_not_rounded_before_the_component(self):
        unit = np.nextafter(0.0, 1.0)
        for first, second, rho, imaginary_units in ((1, 2, 0.25, 1),
                                                     (2, 3, 0.125, 1),
                                                     (1, 2, 1.0, 3)):
            result = math_core.complex_second_order_moments(first * unit, second * unit, rho)
            self.assertEqual(result['pseudo_variance_c'].imag, imaginary_units * unit)
        result = math_core.complex_second_order_moments(unit, 2 * unit, 0.5)
        self.assertEqual(result['pseudo_variance_abs'], 2 * unit)

    def test_moment_magnitude_keeps_cancellation_and_finite_rank_one_boundary(self):
        for first in (1e200, 1e300):
            second = math.nextafter(first, 0.0)
            result = math_core.complex_second_order_moments(first, second)
            self.assertEqual(result['pseudo_variance_abs'], first - second)
        largest = np.finfo(float).max
        for rho in (-1.0, 1.0):
            result = math_core.complex_second_order_moments(largest / 4, largest - largest / 4, rho)
            self.assertEqual(result['pseudo_variance_abs'], largest)

    def test_subnormal_half_unit_is_decided_before_square_root_rounding(self):
        unit = np.nextafter(0.0, 1.0)
        rho = 0.14433756729740646
        # Squaring the exact binary inputs proves this is above half a unit.
        self.assertGreater(48 * Fraction.from_float(rho)**2, 1)
        for sign in (-1, 1):
            result = math_core.complex_second_order_moments(unit, 3 * unit, sign * rho)
            self.assertEqual(result['pseudo_variance_c'].imag, sign * unit)
            with self.assertRaisesRegex(ValueError, 'component'):
                math_core.complex_second_order_moments(unit, unit, sign * 0.25)

    def test_original_exercise_moments_are_retained(self):
        result = math_core.complex_second_order_moments(0.8, 0.2, 0.5)
        self.assertEqual(result['pseudo_variance_abs'], math.sqrt(0.52))


class ReviewedMathematicalContracts(unittest.TestCase):
    def test_column_rotation_products_are_named_without_reversing_application_order(self):
        result = math_core.so3_rotation_matrix_composition_example()
        # Ry first sends ex to -ez; Rz leaves that axis unchanged.
        np.testing.assert_allclose(result['ex_under_r_zyx'], [0, 0, -1], atol=1e-15)
        # Rz first sends ex to ey; Ry leaves that axis unchanged.
        np.testing.assert_allclose(result['ex_under_r_yzx'], [0, 1, 0], atol=1e-15)
        self.assertNotIn('ex_under_rz_then_ry', result)

    def test_real_kalman_rejects_each_complex_covariance_before_losing_its_imaginary_part(self):
        arguments = dict(f_mat=np.eye(2), h_mat=np.eye(2), q_cov=np.eye(2),
                         r_cov=np.eye(2), s_prev=np.zeros(2), p_prev=np.eye(2),
                         z_obs=np.zeros(2))
        hermitian = np.array([[1, .2j], [-.2j, 1]])
        for name in ('q_cov', 'r_cov', 'p_prev'):
            for covariance in (hermitian, np.eye(2, dtype=complex)):
                with self.subTest(name=name, covariance=covariance):
                    with self.assertRaisesRegex(ValueError, 'real-valued'):
                        math_core.kalman_joseph_step(**(arguments | {name: covariance}))

    def test_real_kalman_keeps_a_valid_singular_posterior(self):
        result = math_core.kalman_joseph_step(
            np.eye(2), np.array([[1., 0.]]), np.zeros((2, 2)),
            np.ones((1, 1)), np.zeros(2), np.diag([1., 0.]), np.zeros(1))
        np.testing.assert_allclose(result['posterior_cov_joseph'], np.diag([.5, 0.]))

    def test_nnls_reports_a_nonoptimal_iterate_when_the_budget_is_exhausted(self):
        result = math_core.solve_nnls_projected_gradient(
            np.array([[1., .5], [.5, 1.]]), np.array([1., -.2]), max_iters=1)
        # P^T b=(.9,.3), largest eigenvalue of P^T P=9/4.
        np.testing.assert_allclose(result['q'], [2 / 5, 2 / 15], atol=1e-15)
        self.assertFalse(result['converged'])
        self.assertEqual(result['iterations'], 1)
        self.assertEqual(result['stopping_reason'], 'iteration_limit')
        self.assertAlmostEqual(result['projected_gradient_residual'], 16 / 135)
        self.assertLess(result['dual_min_entry'], 0)

    def test_nnls_normal_case_has_an_independent_kkt_certificate(self):
        result = math_core.solve_nnls_projected_gradient(
            np.array([[1., .5], [.5, 1.]]), np.array([1., -.2]))
        # On the first-axis active face, q1=(P[:,0]^T b)/||P[:,0]||²=.72.
        np.testing.assert_allclose(result['q'], [.72, 0.], atol=1e-14)
        np.testing.assert_allclose(result['dual_multiplier_mu'], [0., .42], atol=1e-14)
        self.assertTrue(result['converged'])
        self.assertEqual(result['stopping_reason'], 'kkt_tolerance')

    def test_nnls_common_problem_scaling_keeps_the_same_scalar_optimum(self):
        for scale in (1., 1e150, 1e155, 1e-155, 1e308, np.nextafter(0., 1.)):
            with self.subTest(scale=scale), warnings.catch_warnings():
                warnings.simplefilter('error', RuntimeWarning)
                result = math_core.solve_nnls_projected_gradient([[scale]], [scale])
            np.testing.assert_array_equal(result['q'], [1.])
            self.assertTrue(result['converged'])
            self.assertEqual(result['objective_value'], 0.)

    def test_nnls_diagnostics_stay_in_original_units_when_representable(self):
        scale = 1e150
        result = math_core.solve_nnls_projected_gradient(
            scale * np.array([[1., .5], [.5, 1.]]), scale * np.array([1., -.2]))
        np.testing.assert_allclose(result['q'], [.72, 0.], atol=1e-14)
        self.assertAlmostEqual(result['objective_value'] / scale / scale, .196)
        self.assertAlmostEqual(result['dual_multiplier_mu'][1] / scale / scale, .42)

    def test_unrepresentable_nnls_problem_or_diagnostics_raise_before_returning_nan(self):
        for matrix, observation in (
                ([[1e-308]], [1e308]), ([[1e308]], [-1e308]),
                (np.diag([1e200, 1.]), [0., 1.]),
                (np.diag([1e200, 1e40]), [0., 1.]),
                ([[1e-170]], [-1e-170])):
            with self.subTest(matrix=matrix), warnings.catch_warnings():
                warnings.simplefilter('error', RuntimeWarning)
                with self.assertRaisesRegex(ValueError, 'floating-point range'):
                    math_core.solve_nnls_projected_gradient(matrix, observation)


if __name__ == '__main__':
    unittest.main()
