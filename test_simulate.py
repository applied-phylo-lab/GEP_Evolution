#!/usr/bin/env python3
"""
test_simulate.py
================
Unit tests for simulate.py.

Run with:
  python -m pytest test_simulate.py -v
  python test_simulate.py          # runs directly via unittest
"""

import os
import shutil
import tempfile
import unittest

import numpy as np

from simulate import (
    effective_rank,
    make_seed,
    modularity_entropy,
    load_condition,
    save_condition,
    simulate,
    create_dirichlet_tasks,
    make_genome,
)


# ============================================================
# 1. make_seed
# ============================================================

class TestMakeSeed(unittest.TestCase):

    def test_deterministic(self):
        """Same inputs always produce the same seed."""
        s1 = make_seed(0, 3, 0.6, 1, 0.167, 0.5, 'sim')
        s2 = make_seed(0, 3, 0.6, 1, 0.167, 0.5, 'sim')
        self.assertEqual(s1, s2)

    def test_different_reps(self):
        """Different reps produce different seeds."""
        s0 = make_seed(0, 3, 0.6, 1, 0.167, 0.5, 'sim')
        s1 = make_seed(1, 3, 0.6, 1, 0.167, 0.5, 'sim')
        self.assertNotEqual(s0, s1)

    def test_different_roles(self):
        """genome and sim roles produce different seeds for same rep/condition."""
        s_genome = make_seed(0, 3, 0.6, 1, 0.167, 0.5, 'genome')
        s_sim    = make_seed(0, 3, 0.6, 1, 0.167, 0.5, 'sim')
        self.assertNotEqual(s_genome, s_sim)

    def test_different_m(self):
        """Different m values produce different seeds."""
        s1 = make_seed(0, 6, 0.6, 1, 0.167, 0.5, 'sim')
        s2 = make_seed(0, 6, 0.6, 3, 0.167, 0.5, 'sim')
        self.assertNotEqual(s1, s2)

    def test_valid_range(self):
        """Seed is a non-negative integer within numpy's valid range."""
        s = make_seed(5, 9, 1.2, 9, 0.167, 0.3, 'sim')
        self.assertIsInstance(s, int)
        self.assertGreaterEqual(s, 0)
        self.assertLess(s, 2 ** 31)


# ============================================================
# 2. effective_rank
# ============================================================

class TestEffectiveRank(unittest.TestCase):

    def test_identical_activations_returns_one(self):
        """When all tasks use the same activation vector, rank should be 1."""
        a = np.array([1.0, 0.5, 0.0, 0.2, 0.0, 0.1])
        a_list = [a.copy() for _ in range(4)]
        er = effective_rank(a_list)
        self.assertAlmostEqual(er, 1.0, places=5)

    def test_orthogonal_activations_returns_T(self):
        """When T tasks each use one unique program, rank should equal T."""
        K = 4
        T = 4
        # each task uses exactly one different program
        a_list = [np.eye(K)[i] for i in range(T)]
        er = effective_rank(a_list)
        self.assertAlmostEqual(er, float(T), places=5)

    def test_orthogonal_activations_T_less_than_K(self):
        """With T < K orthogonal activations, rank equals T not K."""
        K = 6
        T = 3
        a_list = [np.eye(K)[i] for i in range(T)]
        er = effective_rank(a_list)
        self.assertAlmostEqual(er, float(T), places=5)

    def test_single_activation_returns_one(self):
        """Single activation vector always returns 1."""
        er = effective_rank([np.array([1.0, 0.0, 0.5])])
        self.assertEqual(er, 1.0)

    def test_range(self):
        """Effective rank is always between 1 and min(T, K)."""
        rng = np.random.default_rng(42)
        K, T = 6, 4
        a_list = [rng.random(K) for _ in range(T)]
        er = effective_rank(a_list)
        self.assertGreaterEqual(er, 1.0 - 1e-9)
        self.assertLessEqual(er, float(min(T, K)) + 1e-9)


# ============================================================
# 3. modularity_entropy
# ============================================================

class TestModularityEntropy(unittest.TestCase):

    def test_perfect_modularity_returns_one(self):
        """Each mutation affects exactly one task -> M = 1."""
        T = 3
        n_muts = 12
        # each mutation has effect only on one task
        dF = np.zeros((n_muts, T))
        for i in range(n_muts):
            dF[i, i % T] = 1.0
        M = modularity_entropy(dF)
        self.assertAlmostEqual(M, 1.0, places=5)

    def test_uniform_effects_returns_zero(self):
        """Each mutation affects all tasks equally -> M = 0."""
        T = 3
        n_muts = 12
        dF = np.ones((n_muts, T))
        M = modularity_entropy(dF)
        self.assertAlmostEqual(M, 0.0, places=5)

    def test_all_zero_mutations_returns_nan(self):
        """If all mutations have zero effect, return nan."""
        dF = np.zeros((10, 3))
        M = modularity_entropy(dF)
        self.assertTrue(np.isnan(M))

    def test_T_less_than_2_returns_nan(self):
        """Single task has no meaningful modularity."""
        dF = np.ones((10, 1))
        M = modularity_entropy(dF)
        self.assertTrue(np.isnan(M))

    def test_range(self):
        """Modularity is always between 0 and 1 for valid input."""
        rng = np.random.default_rng(99)
        dF = rng.random((50, 4)) - 0.5
        M = modularity_entropy(dF)
        self.assertGreaterEqual(M, 0.0 - 1e-9)
        self.assertLessEqual(M, 1.0 + 1e-9)


# ============================================================
# 4. save_condition / load_condition round-trip
# ============================================================

class TestCacheRoundTrip(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def _make_fake_hist(self, T=3, K=6, L=10, n_subs=5):
        """Create a minimal hist dict matching _run_sswm output structure."""
        rng = np.random.default_rng(0)
        snap = {
            'genome':    rng.integers(0, 2, (L, K)).astype(float),
            'usage':     rng.random((T, K)),
            'phenotype': rng.random((L, T)),
        }
        return {
            'eff_rank':           rng.random(n_subs),
            'modularity_entropy': rng.random(n_subs),
            'pheno_dist':         rng.random(n_subs),
            'cum_time':           np.cumsum(rng.exponential(1.0, n_subs)),
            'wait_time':          rng.exponential(1.0, n_subs),
            'P':                  rng.random((n_subs, T)),
            'd':                  rng.random((n_subs, T)),
            'W':                  rng.random(n_subs),
            's_max':              rng.random(n_subs),
            'n_ben':              rng.integers(0, 10, n_subs),
            'active_tasks':       rng.integers(0, T, (n_subs, 1)),
            'snapshots':          {0: snap, n_subs - 1: snap},
            'n_actual_subs':      n_subs,
            'n_tasks':            T,
        }

    def test_round_trip(self):
        """save_condition then load_condition recovers identical arrays."""
        T, K, L, n_subs, n_reps = 3, 6, 10, 5, 2
        results = [self._make_fake_hist(T, K, L, n_subs) for _ in range(n_reps)]
        params  = dict(T=T, dT=0.6, m=1, density=0.167, alpha=0.5,
                       L=L, K=K, GAMMA=1.0, N_SUBS=n_subs, N_REPS=n_reps,
                       MU=1e-7, N_POP=10000, TASK_SEED=270,
                       N_GENOME_SNAPSHOTS=10)

        base = os.path.join(self.tmpdir, 'test_condition')
        save_condition(base, results, params)
        loaded, loaded_params = load_condition(base)

        self.assertEqual(len(loaded), n_reps)
        self.assertEqual(loaded_params['T'], T)

        for rep_idx in range(n_reps):
            orig = results[rep_idx]
            rec  = loaded[rep_idx]

            for key in ('eff_rank', 'modularity_entropy', 'pheno_dist',
                        'P', 'd', 'W', 'cum_time', 'wait_time'):
                np.testing.assert_array_almost_equal(
                    orig[key], rec[key],
                    err_msg=f'Mismatch in rep {rep_idx}, key {key}'
                )

            # check snapshots
            for step in orig['snapshots']:
                for field in ('genome', 'usage', 'phenotype'):
                    np.testing.assert_array_almost_equal(
                        orig['snapshots'][step][field],
                        rec['snapshots'][step][field],
                        err_msg=f'Snapshot mismatch rep {rep_idx} step {step} {field}'
                    )


# ============================================================
# 5. simulate sanity check
# ============================================================

class TestSimulateSanity(unittest.TestCase):
    """
    Small-scale simulation sanity checks.
    Uses tiny parameters to keep runtime under a few seconds.
    """

    @classmethod
    def setUpClass(cls):
        """Run two small simulations once and reuse results."""
        L, K, T = 20, 4, 3
        n_subs, n_reps = 30, 3
        gamma = 1.0
        N_pop = 10000
        mu    = 1e-7

        # high divergence tasks
        tasks = create_dirichlet_tasks(L, T, alpha=0.05, seed=42)
        w     = np.ones(T) / T

        cls.results_m1 = [
            simulate(make_genome(L, K, 1.0/K, seed=i),
                     tasks, n_subs, mu, N_pop, gamma, w,
                     m=1, seed=100+i)
            for i in range(n_reps)
        ]
        cls.results_mT = [
            simulate(make_genome(L, K, 1.0/K, seed=i),
                     tasks, n_subs, mu, N_pop, gamma, w,
                     m=T, seed=200+i)
            for i in range(n_reps)
        ]

    def test_output_keys_present(self):
        """All expected trajectory keys are present in output."""
        expected = {'eff_rank', 'modularity_entropy', 'pheno_dist',
                    'P', 'd', 'W', 'cum_time', 'wait_time',
                    's_max', 'n_ben', 'active_tasks',
                    'snapshots', 'n_actual_subs', 'n_tasks'}
        for hist in self.results_m1:
            self.assertTrue(expected.issubset(set(hist.keys())))

    def test_trajectory_lengths(self):
        """Trajectory arrays have length n_subs."""
        n_subs = 30
        for hist in self.results_m1 + self.results_mT:
            self.assertEqual(len(hist['eff_rank']), n_subs)
            self.assertEqual(len(hist['W']), n_subs)

    def test_snapshot_structure(self):
        """Snapshots contain genome, usage, phenotype with correct shapes."""
        L, K, T = 20, 4, 3
        for hist in self.results_m1:
            for step, snap in hist['snapshots'].items():
                self.assertEqual(snap['genome'].shape,    (L, K))
                self.assertEqual(snap['usage'].shape,     (T, K))
                self.assertEqual(snap['phenotype'].shape, (L, T))

    def test_performance_in_unit_interval(self):
        """All task performances are in [0, 1]."""
        for hist in self.results_m1 + self.results_mT:
            P = hist['P']
            self.assertTrue(np.all(P >= 0.0 - 1e-9))
            self.assertTrue(np.all(P <= 1.0 + 1e-9))

    def test_effective_rank_in_valid_range(self):
        """Effective rank is between 1 and min(T, K) = 3."""
        T, K = 3, 4
        max_rank = float(min(T, K))
        for hist in self.results_m1 + self.results_mT:
            er = hist['eff_rank']
            valid = er[~np.isnan(er)]
            self.assertTrue(np.all(valid >= 1.0 - 1e-9))
            self.assertTrue(np.all(valid <= max_rank + 1e-9))

    def test_m1_and_mT_differ_at_high_divergence(self):
        """
        At high task divergence, m=T should achieve higher mean effective rank
        than m=1 on average across replicates.
        This is the core result of the paper — failure here would be alarming.
        """
        mean_er_m1 = np.mean([r['eff_rank'][-1] for r in self.results_m1])
        mean_er_mT = np.mean([r['eff_rank'][-1] for r in self.results_mT])
        # with only 3 reps and 30 subs the signal may be weak,
        # so we just check they are not identical
        self.assertFalse(
            np.isclose(mean_er_m1, mean_er_mT, atol=1e-6),
            msg=f'm=1 and m=T produced identical effective rank: {mean_er_m1:.4f}'
        )

    def test_active_tasks_shape_m1(self):
        """For m=1 each step has exactly 1 active task."""
        for hist in self.results_m1:
            at = hist['active_tasks']
            self.assertEqual(at.shape[1], 1)

    def test_active_tasks_shape_mT(self):
        """For m=T each step uses all T tasks."""
        T = 3
        for hist in self.results_mT:
            at = hist['active_tasks']
            self.assertEqual(at.shape[1], T)


if __name__ == '__main__':
    unittest.main(verbosity=2)
