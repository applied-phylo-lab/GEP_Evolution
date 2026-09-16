"""Checks that S1 does not turn incomplete trajectories into plateaus."""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'figures'))
import fig_termination as F


def rep(value=0.2, reason='n_subs', length=2):
    return {'task_dT_realized': 1., 'pheno_dist': np.full(length, value),
            'd': np.full((length, 2), 0.5), 'n_ben': np.ones(length - 1),
            'termination_reason': reason}


class TrajectoryTests(unittest.TestCase):
    def test_three_distinct_regimes(self):
        self.assertEqual([m for m, _, _ in F.regimes(8)], [1, 4, 8])
        self.assertEqual(len({str(ls) for _, ls, _ in F.regimes(8)}), 3)
        self.assertEqual([m for m, _, _ in F.regimes(2)], [1, 2])
        with self.assertRaises(ValueError): F.regimes(7)

    def test_budget_and_safeguard_stops_are_not_extended(self):
        for reason in ['n_subs', 'redraw_cap']:
            _, defined = F.trajectory([rep(reason=reason)], 'differentiation', 5, 400)
            np.testing.assert_array_equal(defined[0], [True, True, False, False, False, False])

    def test_no_beneficial_mutations_endpoint_is_retained(self):
        values, defined = F.trajectory([rep(reason='absorbing')], 'n_ben', 5, 400)
        self.assertTrue(defined.all())
        np.testing.assert_array_equal(values[0], [1/400, 0, 0, 0, 0, 0])

    def test_coverage_threshold(self):
        r = [rep(reason='absorbing') for _ in range(179)] + [rep(reason='redraw_cap') for _ in range(21)]
        mean, se, n = F.summarize_trajectory(r, 'differentiation', 5, 400)
        self.assertTrue(np.isfinite(mean[:2]).all())
        self.assertTrue(np.isnan(mean[2:]).all())
        self.assertTrue(np.isnan(se[2:]).all())
        self.assertEqual(n[-1], 179)

    def test_standard_error_uses_contributing_replicates(self):
        values = np.linspace(0, 1, 180)
        r = [rep(float(v), reason='absorbing') for v in values] + [rep(reason='redraw_cap') for _ in range(20)]
        mean, se, n = F.summarize_trajectory(r, 'differentiation', 5, 400)
        self.assertEqual(n[-1], 180)
        self.assertAlmostEqual(mean[-1], values.mean())
        self.assertAlmostEqual(se[-1], values.std(ddof=1)/np.sqrt(180))

    def test_incomplete_grid_is_rejected(self):
        with self.assertRaises(ValueError): F.validate_grid({}, 8, [.2, .8, 1.4])


if __name__ == '__main__': unittest.main()
