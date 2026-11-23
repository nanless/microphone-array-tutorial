"""Fixed PHD illustration: one shared normalizer, without random-filter replay."""
from copy import deepcopy
import json
import math
from pathlib import Path
import random
import unittest

import numpy as np

from scripts.figure_report_contracts import REPORT_BY_NUMBER
from scripts.portable_figure_report_checks import validate_numeric_report


ROOT = Path(__file__).resolve().parents[1]


class PHDFeasibilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original = json.loads((ROOT / REPORT_BY_NUMBER[21].report_relative).read_bytes())
        cls.grid = np.asarray(cls.original['phd_illustration']['grid_deg'])
        # Independent vectorized densities and nodal trapezoid weights.
        cls.gaussian = (np.exp(-((cls.grid+25)/4)**2/2)/4
                        + np.exp(-((cls.grid-40)/5)**2/2)/5)/math.sqrt(2*math.pi)
        cls.weights = np.concatenate(([(cls.grid[1]-cls.grid[0])/2],
            (cls.grid[2:]-cls.grid[:-2])/2, [(cls.grid[-1]-cls.grid[-2])/2]))
        cls.gaussian_area = float(cls.weights @ cls.gaussian)

    def changed(self, values):
        report = deepcopy(self.original)
        report['phd_illustration']['intensity_targets_per_degree'] = list(map(float,values))
        return report

    def test_current_illustration_remains_a_valid_positive(self):
        validate_numeric_report(21,self.original)

    def test_thirty_six_legal_backgrounds_and_upper_endpoint(self):
        upper = math.nextafter(.0005,0.)
        backgrounds = [np.zeros(400),np.full(400,upper),
            np.asarray([0. if i%2 == 0 else upper for i in range(400)]),
            np.asarray([upper if i == 144 else 0. for i in range(400)]),
            np.asarray([0. if i == 288 else upper for i in range(400)])]
        rng = random.Random(20261123)
        backgrounds += [np.asarray([rng.random()**(j%4+1)*upper for _ in range(400)])
                        for j in range(31)]
        self.assertEqual(len(backgrounds),36)
        # The closed upper envelope also admits its limiting endpoint; it is
        # a bounded-background example, not a claim about a seeded RNG draw.
        backgrounds.append(np.full(400,.0005))
        for index,background in enumerate(backgrounds):
            with self.subTest(background=index):
                raw = self.gaussian+background
                total = float(self.weights @ raw)
                values = 2*raw/total
                self.assertAlmostEqual(float(self.weights @ values),2.,places=14)
                reconstructed = total*values/2-self.gaussian
                self.assertLess(float(np.max(np.abs(reconstructed-background))),1e-16)
                validate_numeric_report(21,self.changed(values))

    def test_separately_legal_peaks_have_no_common_normalizer(self):
        values = 2*self.gaussian/self.gaussian_area
        values[144] += 2*.0005/self.gaussian_area
        values[288] -= 2*.0005/self.gaussian_area
        self.assertAlmostEqual(float(self.weights @ values),2.,places=14)
        # These two points alone impose contradictory bounds on the same T.
        required_lower = 2*self.gaussian[288]/values[288]
        required_upper = 2*(self.gaussian[144]+.0005)/values[144]
        self.assertGreater(required_lower-required_upper,.012)
        with self.assertRaisesRegex(ValueError,'no common normalizer'):
            validate_numeric_report(21,self.changed(values))

    def test_zero_or_unrepresentably_small_tail_cannot_hide_in_absolute_slack(self):
        self.assertGreater(float(np.min(self.gaussian)),1e-59)
        values = 2*self.gaussian/self.gaussian_area
        for invalid in (0.,-0.,math.ulp(0.),1e-308):
            with self.subTest(tail=invalid):
                changed = values.copy()
                changed[0] = invalid
                # Tail changes are too small to alter this floating integral.
                self.assertAlmostEqual(float(self.weights @ changed),2.,places=14)
                with self.assertRaises(ValueError):
                    validate_numeric_report(21,self.changed(changed))

    def test_nonfinite_tail_is_rejected_before_any_normalizer_comparison(self):
        values = 2*self.gaussian/self.gaussian_area
        for invalid in (math.nan,math.inf,-math.inf):
            with self.subTest(tail=invalid):
                changed = values.copy()
                changed[0] = invalid
                with self.assertRaisesRegex(ValueError,'nonfinite'):
                    validate_numeric_report(21,self.changed(changed))


if __name__ == '__main__':
    unittest.main()
