#!/usr/bin/env python3
"""
test_simulate.py
================
Acceptance tests for simulate.py. Standard library plus numpy/scipy only.

    python test_simulate.py

Two tests carry the real risk. The chunking test asserts that changing which
process computes a replicate cannot change results, which rests on every seed
being a pure function of (rep, condition) with no state carried between
replicates. The lazy-enumeration test asserts that computing mutant
performances one task column at a time, on demand, yields bit-identical
trajectories to computing the whole table at every step. Re-run both after any
change to seeding, the worker, or `_MutantPerformance`.

Runtimes are kept short by using small L, K and substitution budgets; the tests
check invariants and equivalences, not scientific behaviour.
"""

import itertools
import sys

import numpy as np

import simulate as S


PASS, FAIL = [], []


def identical(a, b):
    """Exact equality, treating NaN as equal to NaN. `modularity_entropy` is
    NaN off the snapshot schedule, so plain `np.array_equal` would report two
    identical trajectories as different."""
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape:
        return False
    if a.dtype.kind == 'f' and b.dtype.kind == 'f':
        return bool(np.array_equal(a, b, equal_nan=True))
    return bool(np.array_equal(a, b))


def check(name, condition, detail=''):
    (PASS if condition else FAIL).append(name)
    print(f'  {"ok  " if condition else "FAIL"}  {name}'
          + (f'   {detail}' if detail and not condition else ''))


# ------------------------------------------------------------
# Fitness: vectorized rows must agree with the scalar reference
# ------------------------------------------------------------

def test_fitness_rows():
    print('\nfitness_power_mean_rows vs fitness_power_mean')
    rng = np.random.default_rng(0)
    for r in (0.0, 1.0, 2.0, -2.0):
        for trial in range(20):
            T = int(rng.integers(2, 6))
            P = rng.random((15, T))
            if trial % 4 == 0:                     # exercise the P = 0 branch
                P[rng.integers(0, 15), rng.integers(0, T)] = 0.0
            w = np.ones(T) / T
            fast = S.fitness_power_mean_rows(P, w, r)
            slow = np.array([S.fitness_power_mean(row, w, r) for row in P])
            ok = np.allclose(fast, slow, rtol=1e-10, atol=1e-12)
            if not ok:
                check(f'r={r} trial={trial}', False,
                      f'max diff {np.max(np.abs(fast - slow)):.3e}')
                return
    check('all exponents agree, including zero performances', True)


# ------------------------------------------------------------
# Absorbing test: separable form vs brute-force enumeration
# ------------------------------------------------------------

def brute_force_exists(P_wt, P_mut, m, r):
    T = P_wt.size
    w = np.ones(m) / m
    for Sset in itertools.combinations(range(T), m):
        Sset = list(Sset)
        W_wt = S.fitness_power_mean(P_wt[Sset], w, r)
        for row in P_mut:
            if S.selection_coeff(S.fitness_power_mean(row[Sset], w, r), W_wt) > 0:
                return True
    return False


def test_absorbing():
    print('\nbeneficial_subset_exists vs brute-force enumeration')
    rng = np.random.default_rng(1)
    mismatch = 0
    for _ in range(300):
        T = int(rng.integers(2, 7))
        m = int(rng.integers(1, T + 1))
        r = float(rng.choice([0.0, 1.0, 2.0, -2.0]))
        P_wt = rng.uniform(0.05, 0.95, size=T)
        # deltas centred slightly below zero so both outcomes occur often
        P_mut = np.clip(P_wt + rng.normal(-0.004, 0.02, size=(25, T)),
                        1e-6, 1.0)
        if S.beneficial_subset_exists(P_wt, P_mut, m, r) != \
                brute_force_exists(P_wt, P_mut, m, r):
            mismatch += 1
    check('300 random cases agree', mismatch == 0, f'{mismatch} mismatches')

    # m = T must reduce to the plain single-subset check
    P_wt = np.array([0.5, 0.4, 0.3])
    P_mut = np.array([[0.49, 0.39, 0.29], [0.51, 0.41, 0.31]])
    check('m=T reduces to the full-set check',
          S.beneficial_subset_exists(P_wt, P_mut, 3, 0.0) is True)
    check('m=T with only deleterious mutants is absorbing',
          S.beneficial_subset_exists(P_wt, P_mut[:1], 3, 0.0) is False)


# ------------------------------------------------------------
# Chunking equivalence
# ------------------------------------------------------------

def run_reps(rep_range, tasks_stack, T, dT, m, density, alpha, L, K, n_subs):
    out = []
    for rep in rep_range:
        g0 = S.make_genome(L, K, density,
                           seed=S.make_genome_seed(rep, L, K, density))
        h = S.simulate(g0, tasks_stack[rep], n_subs, 1e-7, 10_000, 1.0,
                       np.ones(T) / T, m=m,
                       seed=S.make_sim_seed(rep, T, dT, m, density, alpha),
                       n_genome_snapshots=2)
        h['rep_index'] = rep
        out.append(h)
    return out


def test_chunking():
    print('\nchunked vs unchunked replicate assignment')
    L, K, T, dT, m, density, alpha, n_subs, n_reps = 12, 3, 3, 0.8, 1, 0.25, 1.0, 12, 6
    stack = np.stack([S.create_dirichlet_tasks(L, T, alpha,
                                               S.make_ensemble_seed(i, T, dT))
                      for i in range(n_reps)])

    whole = run_reps(range(n_reps), stack, T, dT, m, density, alpha, L, K, n_subs)
    pieces = []
    for start in range(0, n_reps, 2):
        pieces += run_reps(range(start, start + 2), stack, T, dT, m, density,
                           alpha, L, K, n_subs)
    pieces.sort(key=lambda h: h['rep_index'])

    for key in S.STATE_KEYS + S.EVENT_KEYS:
        same = all(identical(a[key], b[key]) for a, b in zip(whole, pieces))
        check(f'{key} bit-identical', same)

    check('termination reasons identical',
          [h['termination_reason'] for h in whole]
          == [h['termination_reason'] for h in pieces])
    check('replicates differ from one another',
          not identical(whole[0]['P'], whole[1]['P']))


# ------------------------------------------------------------
# Lazy mutant enumeration
# ------------------------------------------------------------

def test_lazy_enumeration():
    print('\nlazy per-column enumeration vs full enumeration')
    L, K, T, gamma = 10, 3, 5, 1.0
    tasks = S.create_dirichlet_tasks(L, T, 1.0, 5)
    genome = S.make_genome(L, K, 0.25, seed=9)

    eager = S._MutantPerformance(genome, tasks, gamma)
    full = eager.full()

    partial = S._MutantPerformance(genome, tasks, gamma)
    partial.ensure([2])
    partial.ensure([0, 2])
    check('columns filled out of order match a full table',
          np.array_equal(partial.full(), full))

    manual = np.empty((L * K, T))
    for r in range(L * K):
        g = genome.copy()
        i, j = divmod(r, K)
        g[i, j] = 1.0 - g[i, j]
        manual[r] = S.compute_performance(g, tasks, gamma)['P']
    check('table matches direct per-mutant computation',
          np.array_equal(full, manual))

    # record_modularity='all' forces a full table at every step; the resulting
    # trajectory must be identical to the lazy default.
    m = 2
    g0 = S.make_genome(L, K, 0.25, seed=4)
    kw = dict(n_subs=25, mu=1e-7, N=10_000, gamma=gamma,
              task_weights=np.ones(T) / T, m=m, seed=17, n_genome_snapshots=2)
    lazy = S.simulate(g0, tasks, record_modularity='snapshots', **kw)
    eagerrun = S.simulate(g0, tasks, record_modularity='all', **kw)

    for key in S.STATE_KEYS + S.EVENT_KEYS:
        if key == 'modularity_entropy':
            continue
        check(f'{key} identical under full enumeration',
              identical(lazy[key], eagerrun[key]))
    check('termination reason identical',
          lazy['termination_reason'] == eagerrun['termination_reason'])
    snap_steps = sorted(lazy['snapshots'])
    check('modularity recorded exactly on the snapshot schedule',
          np.array_equal(np.flatnonzero(np.isfinite(lazy['modularity_entropy'])),
                         np.array(snap_steps)))
    check('recorded modularity values agree',
          np.allclose(lazy['modularity_entropy'][snap_steps],
                      eagerrun['modularity_entropy'][snap_steps]))


def test_pfix_array():
    print('\np_fix_array vs p_fix')
    rng = np.random.default_rng(3)
    s_vals = np.concatenate([rng.uniform(-0.1, 0.5, 200),
                             np.array([0.0, 1e-12, 1e-9, 1e-3, 0.5, 50.0,
                                       np.inf])])
    fast = S.p_fix_array(s_vals, 10_000)
    slow = np.array([S.p_fix(float(x), 10_000) for x in s_vals])
    check('all regimes agree including the neutral and infinite limits',
          np.allclose(fast, slow, rtol=1e-12, atol=1e-15, equal_nan=True))


# ------------------------------------------------------------
# Schema invariants
# ------------------------------------------------------------

def test_schema():
    print('\ntrajectory schema')
    L, K, T, m, n_subs = 12, 3, 4, 2, 15
    tasks = S.create_dirichlet_tasks(L, T, 1.0, 7)
    g0 = S.make_genome(L, K, 0.25, seed=3)
    h = S.simulate(g0, tasks, n_subs, 1e-7, 10_000, 1.0, np.ones(T) / T,
                   m=m, seed=11, n_genome_snapshots=3)

    n_states = h['n_states']
    check('n_subs_realized = n_states - 1', h['n_subs_realized'] == n_states - 1)
    check('state arrays have length n_states',
          all(len(h[k]) == n_states for k in S.STATE_KEYS))
    check('event arrays have length n_states - 1',
          all(len(h[k]) == n_states - 1 for k in S.EVENT_KEYS))
    check('budget is respected', n_states - 1 <= n_subs)
    check('termination reason is recorded',
          h['termination_reason'] in ('n_subs', 'absorbing', 'redraw_cap'))
    check('final state is snapshotted', (n_states - 1) in h['snapshots'])
    check('eff_rank is no longer recorded', 'eff_rank' not in h)
    check('snapshots store the genome only',
          set(h['snapshots'][0]) == {'genome', 'cum_time'})
    check('ep_counts is non-decreasing along states',
          bool(np.all(np.diff(h['ep_counts'], axis=0) >= 0)))

    # every epoch is credited exactly once, successful or not
    total = int(h['ep_counts'][-1].sum())
    accounted = m * (len(h['active_tasks']) + int(h['n_failed_epochs'].sum())
                     + h['n_failed_terminal'])
    check('epoch accounting closes', total == accounted,
          f'{total} != {accounted}')

    # P must match the genome recorded in the final snapshot
    gf = h['snapshots'][n_states - 1]['genome']
    check('final snapshot matches final state',
          np.allclose(S.compute_performance(gf, tasks, 1.0)['P'], h['P'][-1]))

    exp = S.expand_snapshot(h['snapshots'][0], tasks, 1.0)
    check('expand_snapshot recovers the state-0 residuals',
          np.allclose(exp['d'], h['d'][0]))
    check('expand_snapshot returns usage and phenotype',
          exp['usage'].shape == (T, K) and exp['phenotype'].shape == (L, T))


def test_budget_rule():
    print('\nsubstitution budget rule')
    check('E*T/m at T=8, m=1', S.n_subs_rule(8, 1, 100) == 801)
    check('floor applies at m=T', S.n_subs_rule(8, 8, 100) == 401)
    check('T=6, m=1', S.n_subs_rule(6, 1, 100) == 601)
    check('epochs override is honoured',
          S.n_subs_rule(8, 1, 100, epochs_per_task=50, floor=1) == 401)


if __name__ == '__main__':
    test_fitness_rows()
    test_pfix_array()
    test_absorbing()
    test_lazy_enumeration()
    test_chunking()
    test_schema()
    test_budget_rule()
    print(f'\n{len(PASS)} passed, {len(FAIL)} failed')
    if FAIL:
        for name in FAIL:
            print(f'  FAILED: {name}')
    sys.exit(1 if FAIL else 0)