"""Independent mathematical boundary tests for chapter-eight float support."""
from fractions import Fraction
import math
import tempfile
from pathlib import Path
from unittest.mock import patch
import unittest
import numpy as np

from codes.chapters.ch08.core.separation import (
    masked_spatial_covariance, mask_mvdr_2x2, guided_activity_posterior, si_sdr,
)
from codes.chapters.ch08.core.gss_teaching import guided_cacgmm_mvdr
from codes.chapters.ch08.examples import gss_teaching_demo as demo


class NumericalSupportTests(unittest.TestCase):
    def test_scm_subnormal_cross_sum_survives_zero_weight_outlier(self):
        x = np.zeros((1, 2, 7))
        x[0, 0, 0] = x[0, 1, 1] = 1
        x[0, :, 2:6] = 2e-162
        x[0, :, 6] = 1e200
        mask = np.array([[1., 1., 1., 1., 1., 1., 0.]])
        # Raw binary operands are kept in an exact independent outer-product sum.
        expected = float(4*Fraction.from_float(2e-162)**2/6)
        self.assertEqual(expected, np.nextafter(0., 1.))
        actual = masked_spatial_covariance(x, mask)[0]
        self.assertEqual(actual[0, 1], expected)
        self.assertEqual(actual[1, 0], expected)
        np.testing.assert_allclose(actual.diagonal(), [1/6, 1/6], rtol=1e-15)

    def test_scm_large_cancellation_retains_small_cross_residue(self):
        x = np.array([[[1e150, 1e150, 1e-150], [1e150, -1e150, 1e-150]]])
        expected = float(Fraction.from_float(1e-150)**2/3)
        np.testing.assert_allclose(masked_spatial_covariance(x, np.ones((1, 3)))[0, 0, 1],
                                   expected, rtol=5e-16, atol=0)

    def test_scm_final_nonzero_cross_underflow_refused_exact_cancellation_legal(self):
        x = np.array([[[1, 0, 1e-162, 1e200], [0, 1, 1e-162, 1e200]]])
        expected = Fraction.from_float(1e-162)**2/3
        self.assertGreater(expected, 0)
        self.assertEqual(float(expected), 0.)
        with self.assertRaisesRegex(ValueError, 'cross term underflows'):
            masked_spatial_covariance(x, [[1, 1, 1, 0]])
        cancelled = np.array([[[1e150, 1e150, 1e-162, 1e-162],
                               [1e150, -1e150, 1e-162, -1e-162]]])
        self.assertEqual(masked_spatial_covariance(cancelled, np.ones((1, 4)))[0, 0, 1], 0.)

    def test_nonzero_normalization_loss_is_not_zero_statistic(self):
        cases = [
            (np.array([[[1e308, 1e-100, 1e-100], [1e308, 1e-100, -1e-100]]]),
             [[0, .9, .1]], [[0, .1, .9]]),
            (np.array([[[0, 1e200, 1e200], [0, 1e200, -1e200]]]),
             [[1e200, 9e-130, 1e-130]], [[1e200, 1e-130, 9e-130]]),
            (np.array([[[1e200, 1e200, 1e-130], [1e200, -1e200, 0]]]),
             [[.9, .1, 0]], [[.1, .9, 0]]),
        ]
        for x, t, v in cases:
            before = x.copy()
            with self.subTest(x=x), self.assertRaisesRegex(ValueError, 'normalization loses'):
                mask_mvdr_2x2(x, t, v)
            np.testing.assert_array_equal(x, before)

    def test_positive_posterior_underflow_refused_but_disabled_zero_legal(self):
        positive = Fraction.from_float(1e-300)**2
        self.assertGreater(positive, 0)
        self.assertEqual(float(positive/(1+positive)), 0.)
        with self.assertRaisesRegex(ValueError, 'positive posterior underflows'):
            guided_activity_posterior([1e-300, 1], [[1e-300, 1]], [[1]])
        np.testing.assert_array_equal(guided_activity_posterior([1e-300, 1], [[1e-300, 1]], [[0]]), [[0, 1]])
        np.testing.assert_array_equal(guided_activity_posterior([0, 1], [[1, 1]], [[1]]), [[0, 1]])

    def test_principal_direction_and_final_beam_loss_are_explicit_support_errors(self):
        for tiny, message in ((1e-160, 'principal direction'), (1e-150, 'beam sum')):
            x = np.array([[[1e200, 0, 0], [tiny*1e200, 1e200, 1.]]])
            # Rank-one [[1,a],[a,a²]] has EXACT principal direction [1,a].
            # Loaded diag(eps,1+eps) gives a strictly positive second weight
            # a*eps / (1+eps+a²*eps); its final contribution at raw x2=1
            # is mathematically nonzero and itself representable.
            a, eps = Fraction.from_float(tiny), Fraction.from_float(1e-6)
            final = a*eps/(1+eps+a*a*eps)
            self.assertGreater(float(final), 0.)
            before = x.copy()
            with self.assertRaisesRegex(ValueError, message):
                mask_mvdr_2x2(x, [[1, 0, 0]], [[0, 1, 0]])
            np.testing.assert_array_equal(x, before)

    def test_subnormal_roundoff_cancellation_has_explicit_undecidable_state(self):
        for perturb in (False, True):
            scale = 1e-308 if perturb else 1e-309
            value = np.nextafter(np.nextafter(scale, np.inf), np.inf) if perturb else scale
            x = np.array([[[scale, value], [1j*scale, -1j*scale]]])
            y, w, diagnostic = mask_mvdr_2x2(x, [[.9, .1]], [[.1, .9]], return_diagnostics=True)
            positions = diagnostic['rounding_underflow_positions']
            self.assertEqual(diagnostic['rounding_underflow_components'], len(positions))
            if not perturb:
                # Preserve the real runner's missed-diagnostic regression,
                # including an exactly zero normalized contraction.
                self.assertEqual(y[0, 1], 0)
                self.assertEqual(diagnostic['rounding_underflow_components'], 1)
                position = positions[0]
                self.assertEqual((position['frequency_index'], position['frame_index'], position['component']), (0, 1, 'real'))
                self.assertIn('cannot distinguish', position['status'])
            else:
                exact = (Fraction(float(w[0, 0].real))*Fraction(float(value)) -
                         Fraction(float(w[0, 1].imag))*Fraction(float(scale)))
                real_positions = [p for p in positions if
                                  (p['frequency_index'], p['frame_index'], p['component']) == (0, 1, 'real')]
                if y[0, 1].real == 0:
                    # Even when float(exact) is minsub, an explicitly marked
                    # roundoff-scale zero is not a proved mathematical null.
                    self.assertEqual(len(real_positions), 1)
                    self.assertIn('cannot distinguish', real_positions[0]['status'])
                    term_sum = (abs(Fraction(float(w[0, 0].real))*Fraction(float(value))) +
                                abs(Fraction(float(w[0, 1].imag))*Fraction(float(scale))))
                    roundoff_band = (8*Fraction(float(np.finfo(float).eps))*term_sum +
                                     Fraction(float(np.nextafter(0., 1.)))/2)
                    self.assertLessEqual(abs(exact), roundoff_band)
                else:
                    self.assertEqual(y[0, 1].real, float(exact))
                    self.assertLess(abs(y[0, 1].real), np.finfo(float).tiny)
                    self.assertEqual(real_positions, [])
                self.assertEqual(y[0, 1].imag, 0)

    def test_subnormal_cancellation_has_two_exact_rational_backend_controls(self):
        # In the real basis diag(1,j), the actual loaded 2x2 system is
        # [[a,-c],[-c,d]] [r,t]^T=[p,q]^T. Solve it with exact fractions,
        # then permit at most two representational ULPs in its returned pair.
        def backend(shifts):
            def solve(matrix, rhs):
                self.assertEqual(matrix.shape, (2, 2))
                self.assertEqual(matrix[0, 0].imag, 0)
                self.assertEqual(matrix[1, 1].imag, 0)
                self.assertEqual(matrix[0, 1].real, 0)
                self.assertEqual(matrix[1, 0], matrix[0, 1].conjugate())
                self.assertEqual(rhs[0].imag, 0)
                self.assertEqual(rhs[1].real, 0)
                a, c, d, p, q = (Fraction(float(v)) for v in
                                  (matrix[0, 0].real, matrix[0, 1].imag,
                                   matrix[1, 1].real, rhs[0].real, rhs[1].imag))
                det = a*d-c*c
                r, t = float((d*p+c*q)/det), float((c*p+a*q)/det)
                for _ in range(abs(shifts[0])):
                    r = math.nextafter(r, math.inf if shifts[0] > 0 else -math.inf)
                for _ in range(abs(shifts[1])):
                    t = math.nextafter(t, math.inf if shifts[1] > 0 else -math.inf)
                rr, tt = Fraction(r), Fraction(t)
                residual = max(abs(a*rr-c*tt-p), abs(-c*rr+d*tt-q))
                norm_scale = max(abs(a)+abs(c), abs(d)+abs(c))*max(abs(rr), abs(tt)) + max(abs(p), abs(q))
                self.assertLessEqual(residual, 4*Fraction(float(np.finfo(float).eps))*norm_scale)
                return np.array([complex(r, 0), complex(0, t)])
            return solve

        unit = np.nextafter(0., 1.)
        for scale, perturb, shifts in ((1e-309, False, (-1, 0)),
                                       (1e-308, True, (2, -2))):
            value = np.nextafter(np.nextafter(scale, np.inf), np.inf) if perturb else scale
            x = np.array([[[scale, value], [1j*scale, -1j*scale]]])
            with patch.object(np.linalg, 'solve', side_effect=backend(shifts)):
                y, w, diagnostic = mask_mvdr_2x2(x, [[.9, .1]], [[.1, .9]], return_diagnostics=True)
            exact = (Fraction(float(w[0, 0].real))*Fraction(float(value)) -
                     Fraction(float(w[0, 1].imag))*Fraction(float(scale)))
            self.assertNotEqual(exact, 0)
            if not perturb:
                self.assertLess(abs(exact), Fraction(float(unit))/2)
                self.assertEqual(float(exact), 0)
                self.assertEqual(y[0, 1], 0)
                self.assertEqual(diagnostic['rounding_underflow_components'], 1)
                position = diagnostic['rounding_underflow_positions'][0]
                self.assertEqual((position['frequency_index'], position['frame_index'], position['component']), (0, 1, 'real'))
                self.assertIn('cannot distinguish', position['status'])
            else:
                self.assertEqual(float(exact), unit)
                self.assertEqual(y[0, 1], unit)
                self.assertEqual(diagnostic['rounding_underflow_components'], 0)

    def test_gss_tiny_nonzero_row_has_finite_norm_and_remains_valid(self):
        x = np.tile(np.eye(2), (1, 6))[None].astype(complex)
        x[0, :, -1] = 1e-170
        result = guided_cacgmm_mvdr(x, np.ones((12, 2)), iterations=1,
                                    absolute_energy_floor=0, relative_energy_floor=0)
        self.assertTrue(result['valid_points'][0, -1])
        self.assertEqual(result['low_energy_bins'], 0)
        self.assertGreater(result['posterior'][0, 0, -1], 0)

    def test_background_is_summed_without_subtracting_unit_target(self):
        x = np.array([[[1, 1, 0, 0, 1, 0], [0, 0, 1, 1, 0, 1]]], dtype=complex)
        activity = [[1, 0], [1, 0], [0, 1], [0, 1], [1, 1], [1, 1]]
        result = guided_cacgmm_mvdr(x, activity, iterations=1, shape_loading=1e-18)
        # A target-only row has exactly the analytically positive background.
        expected = float(3*Fraction.from_float(1e-18)/(4+Fraction.from_float(1e-18)))
        self.assertAlmostEqual(result['posterior'][0, 2, 0]/expected, 1., places=13)
        self.assertGreater(result['other_scm'][0, 0, 0].real, 0)
        np.testing.assert_allclose(result['other_scm'][0].diagonal(), [expected, 1], rtol=1e-13)

    def test_public_numeric_types_are_not_coerced_from_text_or_bool(self):
        for value in (np.ones((1, 2, 4), dtype=bool), np.full((1, 2, 4), '1')):
            for function in (lambda: masked_spatial_covariance(value, np.ones((1, 4))),
                             lambda: mask_mvdr_2x2(value, np.ones((1, 4)), np.ones((1, 4))),
                             lambda: guided_cacgmm_mvdr(value, np.ones((4, 1)), iterations=1)):
                with self.assertRaises(ValueError):
                    function()
        with self.assertRaises(ValueError):
            guided_cacgmm_mvdr(np.ones((1, 2, 4)), np.full((4, 1), '1'), iterations=1)
        for bad in ('yes', 1, 0, None):
            with self.assertRaises(ValueError):
                si_sdr([1, -1], [1, -1], zero_mean=bad)


class GSSMemberSafetyTests(unittest.TestCase):
    def test_check_flag_requires_boolean_before_model_or_writes(self):
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(demo, 'run_experiment', side_effect=AssertionError('must not run')):
                for bad in ('check', 1, None, np.bool_(True)):
                    with self.assertRaisesRegex(ValueError, 'check must be boolean'):
                        demo.generate(Path(tmp), check=bad)
            self.assertEqual(list(Path(tmp).iterdir()), [])

    def test_real_generate_check_is_readonly_and_real_pcm_mutation_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)/'assets'
            demo.generate(directory)
            before = {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in directory.iterdir()}
            demo.generate(directory, check=True)
            self.assertEqual(before, {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in directory.iterdir()})
            p = directory/'mixture.wav'
            mutated = bytearray(p.read_bytes())
            mutated[444] ^= 1
            p.write_bytes(mutated)
            mutation = {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in directory.iterdir()}
            with self.assertRaisesRegex(ValueError, 'stale'):
                demo.generate(directory, check=True)
            self.assertEqual(mutation, {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in directory.iterdir()})

    def test_extra_missing_and_symlink_preflight_without_running_model(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)/'assets'
            directory.mkdir()
            outside = Path(tmp)/'outside'
            outside.write_bytes(b'untouched')
            link = directory/'source_1.wav'
            link.symlink_to(outside)
            with patch.object(demo, 'run_experiment', side_effect=AssertionError('must not run')):
                for check in (False, True):
                    with self.assertRaisesRegex(ValueError, 'ordinary files'):
                        demo.generate(directory, check=check)
            self.assertEqual(outside.read_bytes(), b'untouched')
            link.unlink()
            (directory/'EXTRA.txt').write_text('extra')
            with patch.object(demo, 'run_experiment', side_effect=AssertionError('must not run')):
                with self.assertRaisesRegex(ValueError, r"asset member set differs: missing=\[.*\], extra=\['EXTRA\.txt'\]"):
                    demo.generate(directory)
                with self.assertRaises(ValueError):
                    demo.generate(directory, check=True)


if __name__ == '__main__':
    unittest.main()
