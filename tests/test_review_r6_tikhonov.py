"""Public ridge-solve range controls with independent scalar truths."""
import math
import unittest
import warnings
from decimal import Decimal, localcontext
import numpy as np
from codes.chapters.ch00.core import chapter00_math as cm


def scalar_ridge(g, y, eps):
    with localcontext() as ctx:
        ctx.prec = 180
        gd, yd, ed = (Decimal.from_float(float(x)) for x in (g, y, eps))
        return float(gd * yd / (gd * gd + ed))


class TikhonovRangeTests(unittest.TestCase):
    def test_large_square_does_not_erase_representable_solution(self):
        g, y = np.eye(2) * 1e160, np.array([3e159, -5e159])
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            out = cm.tikhonov_svd_solve(g, y, 1.)
        np.testing.assert_allclose(out['solution'], [scalar_ridge(1e160, a, 1.) for a in y], rtol=8e-15, atol=0)
        self.assertEqual(caught, [])
        self.assertTrue(math.isnan(out['normal_eq_consistency_error']))

    def test_small_square_contributes_to_minimum_ridge(self):
        g = math.ldexp(1., -540)
        out = cm.tikhonov_svd_solve([[g]], [g], math.ldexp(1., -1074))
        self.assertAlmostEqual(float(out['solution'][0].real), 1/65, places=16)

    def test_complex_unitary_basis_preserves_filtered_modes(self):
        q = np.array([[1, 1j], [1j, 1]]) / math.sqrt(2)
        g = 1e160 * (q @ np.diag([1., .5]) @ q.conj().T)
        y = q @ np.array([3e159, 1e159])
        truth = q @ np.array([scalar_ridge(1e160, 3e159, 1.), scalar_ridge(5e159, 1e159, 1.)])
        np.testing.assert_allclose(cm.tikhonov_svd_solve(g, y, 1.)['solution'], truth, rtol=8e-15, atol=0)

    def test_ridge_rank_and_optional_normal_diagnostic_semantics(self):
        out = cm.tikhonov_svd_solve([[1., 2.]], [1.], 1e-20)
        np.testing.assert_allclose(out['solution'], [.2, .4], rtol=0, atol=1e-14)
        self.assertTrue(math.isinf(out['condition_number_2']))
        self.assertTrue(math.isnan(out['normal_eq_consistency_error']))
        np.testing.assert_allclose(cm.tikhonov_svd_solve([[1., 0], [0, 0]], [3., 4.], .25)['solution'], [2.4, 0], rtol=8e-15, atol=0)
        with self.assertRaises(ValueError):
            cm.tikhonov_svd_solve([[1., 2.]], [1.], 0.)

    def test_unrepresentable_main_solution_is_refused(self):
        with self.assertRaises(ValueError):
            cm.tikhonov_svd_solve([[1e-200]], [1e308], 1e-300)

    def test_subnormal_filter_is_not_rounded_before_large_projection(self):
        u = math.ldexp(1., -1074)
        for g, eps in [(3*u, 2.), (5*u, 3.), (math.ldexp(1., -1022), 2.)]:
            with self.subTest(g=g, eps=eps):
                y = math.ldexp(1., 1023)
                out = cm.tikhonov_svd_solve([[g]], [y], eps)
                np.testing.assert_allclose(out['solution'], [scalar_ridge(g, y, eps)], rtol=8e-15, atol=0)

    def test_subnormal_filter_complex_projection(self):
        u = math.ldexp(1., -1074)
        y = math.ldexp(1., 1021) + 1j * math.ldexp(1., 1022)
        truth = complex(scalar_ridge(3*u, y.real, 2.), scalar_ridge(3*u, y.imag, 2.))
        np.testing.assert_allclose(cm.tikhonov_svd_solve([[3*u]], [y], 2.)['solution'], [truth], rtol=8e-15, atol=0)

    def test_nonrepresentable_factor_explicitly_refused_even_if_product_finite(self):
        # The API returns float64 filter factors; a factor below its range
        # remains an explicit support limit, distinct from a nonzero subnormal.
        with self.assertRaises(ValueError):
            cm.tikhonov_svd_solve([[1e-100]], [1e308], 1e308)
