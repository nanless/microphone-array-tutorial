"""Finite regression controls for the zero-event binomial upper-bound interfaces."""
import math
import unittest
from decimal import Decimal, localcontext

from codes.chapters.ch00.core.chapter00_math import binomial_zero_failure_upper_bound
from codes.chapters.ch00.cross_chapter.tracking_time_exercises import zero_failure_upper_bound


class ZeroFailureUpperBoundRegressionTests(unittest.TestCase):
    @staticmethod
    def decimal_bound(n, tail):
        with localcontext() as context:
            context.prec = 110
            return float(1 - (Decimal.from_float(tail).ln() / Decimal(n)).exp())

    def assert_relative_close(self, actual, expected):
        self.assertGreater(actual, 0.0)
        self.assertTrue(math.isfinite(actual))
        self.assertLessEqual(abs(actual - expected), 8 * math.ulp(expected))

    def test_large_trial_count_does_not_cancel_positive_bound(self):
        expected = 2.995732273553991e-18
        self.assertEqual(expected, self.decimal_bound(10**18, 0.05))
        self.assert_relative_close(binomial_zero_failure_upper_bound(10**18, 0.05), expected)

    def test_near_one_tail_and_two_trials_keep_representable_bound(self):
        alpha = math.nextafter(1.0, 0.0)
        expected = 5.551115123125783e-17
        self.assertEqual(expected, self.decimal_bound(2, alpha))
        self.assert_relative_close(binomial_zero_failure_upper_bound(2, alpha), expected)

    def test_default_tail_acceptance_boundary_matches_decimal(self):
        # Independent Decimal110 fixed-n boundary; no simulation or optional stopping.
        for n, target, passed in [(299573225, 1e-8, False),
                                  (299573226, 1e-8, True),
                                  (2995732272, 1e-9, False),
                                  (2995732273, 1e-9, True)]:
            with self.subTest(n=n, target=target):
                expected = self.decimal_bound(n, 0.05)
                actual = binomial_zero_failure_upper_bound(n, 0.05)
                self.assert_relative_close(actual, expected)
                self.assertEqual(expected <= target, passed)
                self.assertEqual(actual <= target, passed)

    def test_parameterizations_preserve_their_small_inputs(self):
        # Converting alpha to 1-alpha would incorrectly reject the first input.
        self.assertEqual(binomial_zero_failure_upper_bound(1, math.ulp(0.0)), 1.0)
        # Converting confidence to a rounded 1-confidence would erase this input.
        self.assertEqual(zero_failure_upper_bound(1, 1e-300), 1e-300)
        self.assert_relative_close(binomial_zero_failure_upper_bound(100, .05),
                                   self.decimal_bound(100, .05))
        with localcontext() as context:
            context.prec = 110
            confidence = Decimal.from_float(.95)
            expected = float(1 - ((1 - confidence).ln() / Decimal(100)).exp())
        self.assert_relative_close(zero_failure_upper_bound(100, .95), expected)

    def test_validation_and_unsupported_underflow_remain_explicit(self):
        for bad in [True, 0, -1, 2.0]:
            with self.subTest(trials=bad):
                with self.assertRaises(ValueError):
                    binomial_zero_failure_upper_bound(bad, .05)
                with self.assertRaises(ValueError):
                    zero_failure_upper_bound(bad, .95)
        for bad in [True, 0.0, 1.0, complex(.05), float('nan')]:
            with self.subTest(alpha=bad):
                with self.assertRaises(ValueError):
                    binomial_zero_failure_upper_bound(2, bad)
        with self.assertRaises(ValueError):
            binomial_zero_failure_upper_bound(10**400, .05)
        with self.assertRaises(ValueError):
            zero_failure_upper_bound(2, math.ulp(0.0))


if __name__ == '__main__':
    unittest.main()
