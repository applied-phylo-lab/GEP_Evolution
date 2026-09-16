#!/usr/bin/env python3
"""
simulate.py
===========
SSWM simulation of gene expression program evolution under variable
simultaneity of selection.

A genome is a binary matrix G in {0,1}^(L x K) mapping L loci onto K programs.
For each of T tasks it deploys the non-negative program combination that best
approximates that task's optimum, and performance decreases with the residual.
Each selective epoch draws m of the T tasks and fitness is a power mean over
them: m = 1 is sequential selection, m = T simultaneous.

Fixation is haploid: p_fix(s) = (1 - e^-2s) / (1 - e^-2Ns), neutral limit 1/N,
mutation supply N*mu.

`termination_reason`, per replicate:
  'n_subs'      the substitution budget was exhausted
  'absorbing'   no single-bit mutant is beneficial under any admissible active
                subset
  'redraw_cap'  MAX_REDRAWS consecutive draws produced nothing

When m < T an empty beneficial set under the drawn subset is not termination:
the epoch is consumed, the subset redrawn, and the walk continues. Failed
epochs still count toward each task's selective exposure.

Substitution budget: S(T, m) = max(floor, round(E*T/m) + 1), E defaulting to L.

Seeding: replicate i has its own initial genome, seeded on (i, L, K, density),
and its own task ensemble, seeded on (i, T, dT). Neither depends on m.

Array lengths: STATE arrays have n_subs_realized + 1 entries, index k being the
genotype after k substitutions. EVENT arrays have one fewer, index k describing
the substitution from state k to state k+1. Arrays are not forward-filled.

Cache layout:
  tasks_L{L}_v{V}/taskens_T{T}.npz
  L{L}_K{K}_gamma{G}_fr{R}_v{V}/density{D}/sim_T{T}_dT{dT}_m{m}_alpha{A}.npz

Usage:
  python simulate.py --dry_run
  python simulate.py --T 2 4 8 --m 1 2 4 --n_reps 100
"""

import argparse
import hashlib
import json
import os
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import nnls


SCHEMA_VERSION = 2

DEFAULTS = dict(
    N_POP=10_000,
    MU=1e-7,
    L=100,
    K=4,
    GAMMA=1.0,
    FITNESS_R=0.0,
    N_REPS=200,
    T_VALUES=[2, 4, 6, 8],
    M_VALUES=None,
    TASK_DIVERGENCES=[0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4],
    GENOME_DENSITIES=[0.25],
    TASK_SEED=270,
    CACHE_DIR='simulation_cache',
    N_WORKERS=None,              # None -> os.cpu_count()
    N_GENOME_SNAPSHOTS=10,
    TASK_SAMPLING='random',      # 'random' | 'sliding'
    N_SUBS=None,                 # None -> use the rule below
    N_SUBS_EPOCHS=None,          # target epochs per task; None -> L
    N_SUBS_FLOOR=401,
    MAX_REDRAWS=100,
    CHUNK_SIZE=10,               # replicates per work item
    RECORD_MODULARITY='snapshots',   # 'snapshots' | 'all' | 'none'
)


# ============================================================
# 1. SEEDS
# ============================================================

def _seed_from(key: str) -> int:
    return int(hashlib.md5(key.encode()).hexdigest(), 16) % (2 ** 31)


def make_genome_seed(rep: int, L: int, K: int, density: float) -> int:
    """Initial genome seed; independent of T, dT and m."""
    return _seed_from(f'genome_{rep}_{L}_{K}_{density:.10g}')


def make_ensemble_seed(rep: int, T: int, dT: float) -> int:
    """Task ensemble seed; independent of m."""
    return _seed_from(f'ensemble_{rep}_{T}_{dT:.10g}')


def make_sim_seed(rep: int, T: int, dT: float, m: int,
                  density: float, alpha: float) -> int:
    """Stochastic path seed: task draws, mutation choice, waiting times."""
    return _seed_from(f'sim_{rep}_{T}_{dT:.10g}_{m}_{density:.10g}_{alpha:.10g}')


# ============================================================
# 2. TASK ENSEMBLES
# ============================================================

def create_dirichlet_tasks(L: int, T: int, alpha: float, seed: int) -> np.ndarray:
    """(L, T) array of unit-norm non-negative task optima."""
    rng = np.random.default_rng(seed)
    tasks = rng.dirichlet(np.ones(L) * alpha, size=T).T
    for i in range(T):
        tasks[:, i] /= np.linalg.norm(tasks[:, i])
    return tasks


def compute_mean_pairwise_distance(tasks: np.ndarray) -> float:
    T = tasks.shape[1]
    dists = [np.linalg.norm(tasks[:, i] - tasks[:, j])
             for i, j in combinations(range(T), 2)]
    return float(np.mean(dists)) if dists else 0.0


def find_alpha_for_target(target_dist: float, L: int, T: int,
                          n_seeds: int = 30, tol: float = 0.02) -> float:
    """Bisect in log-space for the alpha whose expected mean pairwise task
        divergence matches `target_dist`.
    """
    def mean_dist(alpha):
        return np.mean([compute_mean_pairwise_distance(
            create_dirichlet_tasks(L, T, alpha, s)) for s in range(n_seeds)])

    lo, hi = 0.01, 100.0
    for _ in range(50):
        mid = np.sqrt(lo * hi)
        val = mean_dist(mid)
        if abs(val - target_dist) < tol:
            return mid
        if val > target_dist:
            lo = mid
        else:
            hi = mid
    return mid


def build_task_ensembles(T: int, L: int, task_divergences: List[float],
                         n_ensembles: int) -> Tuple[Dict, Dict, Dict]:
    """Calibrate one alpha per dT and draw `n_ensembles` ensembles at each.

        Returns (alpha_map, ensembles, realized); ensembles[dT] is
        (n_ensembles, L, T) and realized[dT] holds each ensemble's mean pairwise
        divergence.
    """
    alpha_map, ensembles, realized = {}, {}, {}

    for dT in task_divergences:
        alpha = find_alpha_for_target(dT, L, T)
        stack = np.stack([
            create_dirichlet_tasks(L, T, alpha, make_ensemble_seed(rep, T, dT))
            for rep in range(n_ensembles)
        ])
        alpha_map[dT] = alpha
        ensembles[dT] = stack
        realized[dT] = np.array([compute_mean_pairwise_distance(stack[i])
                                 for i in range(n_ensembles)])

    return alpha_map, ensembles, realized


# ============================================================
# 3. GENOME
# ============================================================

def make_genome(L: int, K: int, density: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return (rng.random((L, K)) < density).astype(float)


# ============================================================
# 4. PERFORMANCE
# ============================================================

def _task_performance(genome: np.ndarray, task: np.ndarray, gamma: float):
    """(d, P, a, z) for one task: residual distance of the optimal non-negative
        deployment, performance, activation vector, expressed phenotype.
    """
    a, _ = nnls(genome, task)
    z = genome @ a
    d = float(np.linalg.norm(task - z))
    return d, max(0.0, 1.0 - d) ** gamma, a, z


def compute_performance(genome: np.ndarray, tasks: np.ndarray,
                        gamma: float) -> Dict:
    """Optimal non-negative deployment of `genome` against each task.

        Returns d, P, a, usage (T, K) and phenotype (L, T).
    """
    n_tasks = tasks.shape[1]
    d = np.zeros(n_tasks)
    P = np.zeros(n_tasks)
    a_list = []
    phenotypes = np.zeros((genome.shape[0], n_tasks))

    for e in range(n_tasks):
        d[e], P[e], a, z = _task_performance(genome, tasks[:, e], gamma)
        a_list.append(a)
        phenotypes[:, e] = z

    return {'d': d, 'P': P, 'a': a_list,
            'usage': np.array(a_list), 'phenotype': phenotypes}


class _MutantPerformance:
    """Lazily filled (L*K, T) table of single-bit mutant performances.

        Row i mutates entry divmod(i, K); columns are computed on demand.
    """

    __slots__ = ('_genome', '_tasks', '_gamma', '_P', '_have', '_L', '_K', '_T')

    def __init__(self, genome: np.ndarray, tasks: np.ndarray, gamma: float):
        self._genome = genome
        self._tasks = tasks
        self._gamma = gamma
        self._L, self._K = genome.shape
        self._T = tasks.shape[1]
        self._P = np.empty((self._L * self._K, self._T), dtype=float)
        self._have = np.zeros(self._T, dtype=bool)

    def ensure(self, cols: Sequence[int]) -> None:
        todo = [int(j) for j in np.unique(np.asarray(cols)) if not self._have[j]]
        if not todo:
            return
        row = 0
        for i in range(self._L):
            for k in range(self._K):
                mut = self._genome.copy()
                mut[i, k] = 1.0 - mut[i, k]
                for j in todo:
                    _, P, _, _ = _task_performance(mut, self._tasks[:, j],
                                                   self._gamma)
                    self._P[row, j] = P
                row += 1
        for j in todo:
            self._have[j] = True

    def block(self, cols: Sequence[int]) -> np.ndarray:
        """(L*K, len(cols)) copy of the requested columns, computed as needed."""
        self.ensure(cols)
        return self._P[:, np.asarray(cols)]

    def full(self) -> np.ndarray:
        """(L*K, T) table, computing any columns not yet filled."""
        self.ensure(range(self._T))
        return self._P


def make_snapshot(genome: np.ndarray, cum_time: float) -> Dict:
    """Genome and cumulative time; see `expand_snapshot`."""
    return {'genome': genome.copy(), 'cum_time': float(cum_time)}


def expand_snapshot(snapshot: Dict, tasks: np.ndarray, gamma: float) -> Dict:
    """Recover program usage, expressed phenotypes, residuals and performances
    from a stored snapshot."""
    info = compute_performance(np.asarray(snapshot['genome'], dtype=float),
                               tasks, gamma)
    return {'genome': snapshot['genome'], 'cum_time': snapshot['cum_time'],
            'usage': info['usage'], 'phenotype': info['phenotype'],
            'd': info['d'], 'P': info['P']}


# ============================================================
# 5. FITNESS
# ============================================================

def fitness_power_mean(P: np.ndarray, w: np.ndarray, r: float) -> float:
    """Weighted power mean of task performances with exponent r.

        r = 0 is the geometric mean. For r <= 0 a non-positive performance gives 0.
    """
    if r == 0.0:
        if np.any(P <= 0):
            return 0.0
        return float(np.exp(np.sum(w * np.log(P))))
    if np.any(P <= 0) and r < 0:
        return 0.0
    return float(np.power(np.sum(w * np.power(P, r)), 1.0 / r))


def fitness_power_mean_rows(P_rows: np.ndarray, w: np.ndarray,
                            r: float) -> np.ndarray:
    """Vectorized `fitness_power_mean` over the rows of P_rows."""
    P_rows = np.asarray(P_rows, dtype=float)
    if r == 0.0:
        bad = np.any(P_rows <= 0, axis=1)
        safe = np.where(P_rows <= 0, 1.0, P_rows)
        out = np.exp(np.sum(w * np.log(safe), axis=1))
        return np.where(bad, 0.0, out)
    if r < 0:
        bad = np.any(P_rows <= 0, axis=1)
        safe = np.where(P_rows <= 0, 1.0, P_rows)
        out = np.power(np.sum(w * np.power(safe, r), axis=1), 1.0 / r)
        return np.where(bad, 0.0, out)
    return np.power(np.sum(w * np.power(P_rows, r), axis=1), 1.0 / r)


# ============================================================
# 6. POPULATION GENETICS  (haploid; see module docstring)
# ============================================================

def selection_coeff(W_mut: float, W_wt: float) -> float:
    if W_wt <= 0:
        return float('inf') if W_mut > 0 else 0.0
    return (W_mut - W_wt) / W_wt


def selection_coeff_array(W_mut: np.ndarray, W_wt: float) -> np.ndarray:
    """Vectorized `selection_coeff`. A non-positive wild-type fitness makes any
        positive mutant infinitely favoured.
    """
    W_mut = np.asarray(W_mut, dtype=float)
    if W_wt <= 0:
        return np.where(W_mut > 0, np.inf, 0.0)
    return (W_mut - W_wt) / W_wt


def p_fix(s: float, N: int) -> float:
    """Haploid fixation probability, (1 - e^-2s) / (1 - e^-2Ns); neutral limit 1/N."""
    if s <= 0:
        return 0.0
    Ns = N * s
    if Ns > 100:
        return 1.0 - np.exp(-2.0 * s)
    if Ns < 1e-6:
        return (1.0 / N) * (1.0 + (N - 1) * s)
    return (1.0 - np.exp(-2.0 * s)) / (1.0 - np.exp(-2.0 * Ns))


def p_fix_array(s: np.ndarray, N: int) -> np.ndarray:
    """Vectorized `p_fix`."""
    s = np.asarray(s, dtype=float)
    out = np.zeros_like(s)
    pos = s > 0
    if not np.any(pos):
        return out

    sp = s[pos]
    Ns = N * sp
    res = np.empty_like(sp)

    large = Ns > 100
    small = Ns < 1e-6
    mid = ~(large | small)

    res[large] = 1.0 - np.exp(-2.0 * sp[large])
    res[small] = (1.0 / N) * (1.0 + (N - 1) * sp[small])
    if np.any(mid):
        res[mid] = ((1.0 - np.exp(-2.0 * sp[mid]))
                    / (1.0 - np.exp(-2.0 * Ns[mid])))

    out[pos] = res
    return out


# ============================================================
# 7. ABSORBING-STATE TEST
# ============================================================

def beneficial_subset_exists(P_wt: np.ndarray, P_mut: np.ndarray,
                             m: int, r: float) -> bool:
    """Does any single-bit mutant beat the wild type under at least one
        admissible active subset of size m? Exact, not sampled.

        Assumes uniform task weights. The power mean is then separable over the
        active set, so the best subset is the m largest (r >= 0) or m smallest
        (r < 0) per-task deltas, at O(L*K*T log T) instead of C(T, m) subsets. When
        r <= 0 and a performance is non-positive the separable form fails and all
        C(T, m) subsets are enumerated instead. `P_mut` is the full (L*K, T) table.
    """
    T = P_wt.size
    m = min(m, T)

    degenerate = (r <= 0) and (P_wt.min() <= 0 or P_mut.min() <= 0)

    if not degenerate:
        if r == 0.0:
            delta = np.log(P_mut) - np.log(P_wt)
        else:
            delta = np.power(P_mut, r) - np.power(P_wt, r)

        if r < 0:
            best = np.partition(delta, m - 1, axis=1)[:, :m].sum(axis=1)
            return bool(np.any(best < 0.0))
        best = np.partition(delta, T - m, axis=1)[:, T - m:].sum(axis=1)
        return bool(np.any(best > 0.0))

    w = np.ones(m) / m
    for S in combinations(range(T), m):
        S = list(S)
        W_wt = fitness_power_mean(P_wt[S], w, r)
        W_mut = fitness_power_mean_rows(P_mut[:, S], w, r)
        if np.any(selection_coeff_array(W_mut, W_wt) > 0):
            return True
    return False


# ============================================================
# 8. METRICS
# ============================================================

def modularity_entropy(dF: np.ndarray) -> float:
    """Mutational modularity, 1 - mean(H) / log(T), over each mutation's absolute
        effects across tasks. dF is (L*K, T); mutations with no total effect are
        excluded.
    """
    T = dF.shape[1]
    if T < 2:
        return float('nan')
    abs_dF = np.abs(dF)
    row_sums = abs_dF.sum(axis=1, keepdims=True)
    valid = row_sums.flatten() > 1e-12
    if not np.any(valid):
        return float('nan')
    p = abs_dF[valid] / row_sums[valid]
    with np.errstate(divide='ignore', invalid='ignore'):
        log_p = np.where(p > 0, np.log(p), 0.0)
    H = -np.sum(p * log_p, axis=1)
    H_max = np.log(T)
    return float(1.0 - np.mean(H) / H_max) if H_max > 0 else float('nan')


def mean_pairwise_phenotype_distance(phenotype: np.ndarray) -> float:
    """Mean pairwise Euclidean distance between the columns of an (L, T) array."""
    T = phenotype.shape[1]
    if T < 2:
        return 0.0
    return float(np.mean([np.linalg.norm(phenotype[:, i] - phenotype[:, j])
                          for i, j in combinations(range(T), 2)]))


# ============================================================
# 9. SSWM ENGINE
# ============================================================

STATE_KEYS = ('pheno_dist', 'cum_time', 'modularity_entropy',
              'P', 'd', 'ep_counts')
EVENT_KEYS = ('W', 'wait_time', 's_max', 'n_ben', 'n_failed_epochs',
              'active_tasks')


def _draw_active(rng, T: int, m: int, epoch: int, task_sampling: str):
    if m == T:
        return np.arange(T)
    if task_sampling == 'sliding':
        start = epoch % T
        return np.array([i % T for i in range(start, start + m)])
    return rng.choice(T, size=m, replace=False)


def _run_sswm(genome_init: np.ndarray, tasks: np.ndarray, n_subs: int,
              mu: float, N: int, gamma: float, m: int, w: np.ndarray,
              seed: int, n_genome_snapshots: int = 10,
              task_sampling: str = 'random', fitness_r: float = 0.0,
              max_redraws: int = 100,
              record_modularity: str = 'snapshots') -> Dict:
    """Evolve one replicate. See the module docstring for termination reasons and
        the state/event array convention.
    """
    rng = np.random.default_rng(seed)
    L, K = genome_init.shape
    T = tasks.shape[1]

    state = {k: [] for k in STATE_KEYS}
    event = {k: [] for k in EVENT_KEYS}

    ep = np.zeros(T, dtype=np.int64)      # cumulative selective epochs per task
    snapshots: Dict[int, Dict] = {}
    snap_interval = max(n_subs, 1) / (n_genome_snapshots + 1)
    next_snap = 0.0

    genome = genome_init.copy()
    cumtime = 0.0
    epoch_index = 0                       # advances on failed epochs too, so a
                                          # sliding window cannot stall in place
    termination = None
    n_failed_terminal = 0

    for step in range(n_subs + 1):
        info = compute_performance(genome, tasks, gamma)
        P_wt = info['P']
        mutants = _MutantPerformance(genome, tasks, gamma)

        is_snap = step >= next_snap
        if is_snap:
            snapshots[step] = make_snapshot(genome, cumtime)
            next_snap += snap_interval

        if record_modularity == 'all' or (record_modularity == 'snapshots'
                                          and is_snap):
            mod = modularity_entropy(mutants.full() - P_wt)
        else:
            mod = np.nan

        state['pheno_dist'].append(
            mean_pairwise_phenotype_distance(info['phenotype']))
        state['cum_time'].append(cumtime)
        state['modularity_entropy'].append(mod)
        state['P'].append(P_wt.copy())
        state['d'].append(info['d'].copy())
        state['ep_counts'].append(ep.copy())

        if step == n_subs:
            termination = 'n_subs'
            break

        # --- attempt one substitution, redrawing the active subset as needed
        drawn = None
        n_failed = 0
        checked_absorbing = False

        while True:
            active = _draw_active(rng, T, m, epoch_index, task_sampling)
            epoch_index += 1
            ep[active] += 1

            w_act = w[active] / w[active].sum()
            W_wt = fitness_power_mean(P_wt[active], w_act, fitness_r)
            W_mut = fitness_power_mean_rows(mutants.block(active), w_act,
                                            fitness_r)
            s_vals = selection_coeff_array(W_mut, W_wt)
            ben = np.flatnonzero(s_vals > 0)

            if ben.size > 0:
                drawn = (active, W_wt, s_vals, ben)
                break

            n_failed += 1
            if not checked_absorbing:
                if not beneficial_subset_exists(P_wt, mutants.full(), m,
                                                fitness_r):
                    termination = 'absorbing'
                    break
                checked_absorbing = True
            if n_failed >= max_redraws:
                termination = 'redraw_cap'
                break

        if drawn is None:
            n_failed_terminal = n_failed
            break

        active, W_wt, s_vals, ben = drawn
        pfix = p_fix_array(s_vals[ben], N)
        lambda_tot = N * mu * pfix.sum()
        wait = rng.exponential(1.0 / lambda_tot)
        cumtime += wait

        chosen = int(ben[rng.choice(len(ben), p=pfix / pfix.sum())])
        i, j = divmod(chosen, K)

        event['active_tasks'].append(active.copy())
        event['W'].append(W_wt)
        event['wait_time'].append(wait)
        event['s_max'].append(float(s_vals[ben].max()))
        event['n_ben'].append(int(ben.size))
        event['n_failed_epochs'].append(n_failed)

        genome = genome.copy()
        genome[i, j] = 1.0 - genome[i, j]

    n_states = len(state['P'])
    final_step = n_states - 1
    if final_step not in snapshots:
        snapshots[final_step] = make_snapshot(genome, cumtime)

    hist = {k: np.array(v) for k, v in state.items()}
    hist.update({k: np.array(v) for k, v in event.items()})
    hist['snapshots'] = snapshots
    hist['n_states'] = n_states
    hist['n_subs_realized'] = n_states - 1
    hist['n_tasks'] = T
    hist['termination_reason'] = termination
    hist['n_failed_terminal'] = n_failed_terminal
    hist['task_dT_realized'] = compute_mean_pairwise_distance(tasks)
    return hist


def simulate(genome_init: np.ndarray, tasks: np.ndarray, n_subs: int,
             mu: float, N: int, gamma: float, task_weights: np.ndarray,
             m: int, seed: int, n_genome_snapshots: int = 10,
             task_sampling: str = 'random', fitness_r: float = 0.0,
             max_redraws: int = 100,
             record_modularity: str = 'snapshots') -> Dict:
    """Public entry point. `task_weights` is normalized internally."""
    w = task_weights / task_weights.sum()
    return _run_sswm(genome_init, tasks, n_subs, mu, N, gamma, m=m, w=w,
                     seed=seed, n_genome_snapshots=n_genome_snapshots,
                     task_sampling=task_sampling, fitness_r=fitness_r,
                     max_redraws=max_redraws,
                     record_modularity=record_modularity)


# ============================================================
# 10. SUBSTITUTION BUDGET
# ============================================================

def n_subs_rule(T: int, m: int, L: int, epochs_per_task: Optional[int] = None,
                floor: int = 401) -> int:
    """Substitution budget for one condition; `epochs_per_task` defaults to L."""
    E = L if epochs_per_task is None else epochs_per_task
    return max(floor, int(round(E * T / m)) + 1)


# ============================================================
# 11. CACHE
# ============================================================

def _param_root(cache_dir: str, L: int, K: int, gamma: float,
                fitness_r: float) -> str:
    folder = os.path.join(
        cache_dir, f'L{L}_K{K}_gamma{gamma}_fr{fitness_r}_v{SCHEMA_VERSION}')
    os.makedirs(folder, exist_ok=True)
    return folder


def _density_root(cache_dir: str, L: int, K: int, gamma: float,
                  fitness_r: float, density: float) -> str:
    folder = os.path.join(_param_root(cache_dir, L, K, gamma, fitness_r),
                          f'density{density:.4f}')
    os.makedirs(folder, exist_ok=True)
    return folder


def task_cache_path(cache_dir: str, L: int, T: int) -> str:
    """Path to the task ensembles for (L, T), shared across K, gamma and
        fitness_r.
    """
    folder = os.path.join(cache_dir, f'tasks_L{L}_v{SCHEMA_VERSION}')
    os.makedirs(folder, exist_ok=True)
    return os.path.join(folder, f'taskens_T{T}')


def sim_cache_path(cache_dir: str, L: int, K: int, gamma: float,
                   fitness_r: float, density: float, T: int, dT: float,
                   m: int, alpha: float) -> str:
    return os.path.join(
        _density_root(cache_dir, L, K, gamma, fitness_r, density),
        f'sim_T{T}_dT{dT}_m{m}_alpha{alpha:.4f}')


def save_task_ensembles(base_path: str, ensembles: Dict, alpha_map: Dict,
                        realized: Dict, T: int, L: int, n_ensembles: int):
    arrays = {f'tasks_dT{dT}': arr for dT, arr in ensembles.items()}
    arrays.update({f'realized_dT{dT}': r for dT, r in realized.items()})
    np.savez_compressed(base_path + '.npz', **arrays)
    with open(base_path + '_meta.json', 'w') as f:
        json.dump({'T': T, 'L': L, 'n_ensembles': n_ensembles,
                   'schema_version': SCHEMA_VERSION,
                   'dT_values': list(ensembles.keys()),
                   'alpha_map': {str(k): v for k, v in alpha_map.items()},
                   'realized_dT_mean': {str(k): float(np.mean(v))
                                        for k, v in realized.items()},
                   'realized_dT_sd': {str(k): float(np.std(v, ddof=1))
                                      for k, v in realized.items()}},
                  f, indent=2)


def load_task_ensembles(base_path: str) -> Tuple[Dict, Dict, Dict]:
    """Returns (ensembles, realized, meta). ensembles[dT] is (n_ensembles, L, T)."""
    data = np.load(base_path + '.npz')
    with open(base_path + '_meta.json') as f:
        meta = json.load(f)
    ensembles = {dT: data[f'tasks_dT{dT}'] for dT in meta['dT_values']}
    realized = {dT: data[f'realized_dT{dT}'] for dT in meta['dT_values']}
    return ensembles, realized, meta


def _flatten_snapshots(snapshots: Dict[int, Dict]) -> Dict[str, np.ndarray]:
    """Genomes are stored as uint8; usage and phenotype are recomputed by
        `expand_snapshot`.
    """
    flat = {}
    for step, snap in snapshots.items():
        flat[f'snap_{step}_genome'] = snap['genome'].astype(np.uint8)
        flat[f'snap_{step}_cum_time'] = np.array([snap['cum_time']])
    return flat


def _unflatten_snapshots(arrays: Dict, step_keys: List[int]) -> Dict[int, Dict]:
    return {step: {'genome': arrays[f'snap_{step}_genome'].astype(float),
                   'cum_time': float(arrays[f'snap_{step}_cum_time'][0])}
            for step in step_keys}


def save_condition(base_path: str, results: List[Dict], params: Dict):
    """Write all replicates for one condition: trajectory arrays and flattened
        snapshots to npz, per-replicate scalars to json.
    """
    arrays, rep_meta = {}, []

    for slot, hist in enumerate(results):
        prefix = f'rep{slot}_'
        for key in STATE_KEYS + EVENT_KEYS:
            arrays[prefix + key] = hist[key]
        for k, v in _flatten_snapshots(hist['snapshots']).items():
            arrays[prefix + k] = v
        rep_meta.append({
            'rep_index': int(hist['rep_index']),
            'n_states': int(hist['n_states']),
            'n_subs_realized': int(hist['n_subs_realized']),
            'n_tasks': int(hist['n_tasks']),
            'termination_reason': hist['termination_reason'],
            'n_failed_terminal': int(hist['n_failed_terminal']),
            'task_dT_realized': float(hist['task_dT_realized']),
            'snapshot_steps': sorted(hist['snapshots'].keys()),
        })

    np.savez_compressed(base_path + '.npz', **arrays)
    with open(base_path + '_meta.json', 'w') as f:
        json.dump({'params': params, 'n_reps': len(results),
                   'schema_version': SCHEMA_VERSION, 'rep_meta': rep_meta},
                  f, indent=2)


def load_condition(base_path: str,
                   keep: Optional[Sequence[str]] = None
                   ) -> Tuple[List[Dict], Dict]:
    """Load a saved condition as (results, params), sorted by rep_index.

        `keep` restricts which trajectory arrays are read; None reads everything,
        including snapshots.
    """
    npz = np.load(base_path + '.npz', allow_pickle=False)
    with open(base_path + '_meta.json') as f:
        meta = json.load(f)

    if meta.get('schema_version') != SCHEMA_VERSION:
        raise ValueError(
            f'{base_path}: schema version {meta.get("schema_version")} '
            f'!= {SCHEMA_VERSION}')

    wanted = (STATE_KEYS + EVENT_KEYS if keep is None
              else tuple(k for k in STATE_KEYS + EVENT_KEYS if k in keep))
    want_snaps = keep is None or 'snapshots' in keep

    results = []
    for slot in range(meta['n_reps']):
        prefix = f'rep{slot}_'
        info = meta['rep_meta'][slot]
        hist = {key: npz[prefix + key] for key in wanted}
        if want_snaps:
            snap_arrays = {k[len(prefix):]: npz[k] for k in npz.files
                           if k.startswith(prefix + 'snap_')}
            hist['snapshots'] = _unflatten_snapshots(snap_arrays,
                                                     info['snapshot_steps'])
        hist.update({k: info[k] for k in
                     ('rep_index', 'n_states', 'n_subs_realized', 'n_tasks',
                      'termination_reason', 'n_failed_terminal',
                      'task_dT_realized')})
        results.append(hist)

    results.sort(key=lambda h: h['rep_index'])
    return results, meta['params']


def cached_condition_spec(base_path: str) -> Optional[Tuple[int, int]]:
    """(N_SUBS, n_reps) of a cached condition, or None if unreadable."""
    try:
        with open(base_path + '_meta.json') as f:
            meta = json.load(f)
        if meta.get('schema_version') != SCHEMA_VERSION:
            return None
        return int(meta['params']['N_SUBS']), int(meta['n_reps'])
    except (OSError, KeyError, ValueError):
        return None


# ============================================================
# 12. WORKER
# ============================================================

def _worker(args: tuple):
    """Run one chunk of replicates for one condition. Every seed is a pure
        function of (rep, condition), so chunking does not affect results.
    """
    (T, dT, m, density, alpha, tasks_chunk, rep_start, n_subs, mu, N_pop,
     gamma, K, L, task_weights, n_snaps, task_sampling, fitness_r,
     max_redraws, record_modularity) = args

    results = []
    for offset in range(tasks_chunk.shape[0]):
        rep = rep_start + offset
        g0 = make_genome(L, K, density, seed=make_genome_seed(rep, L, K, density))
        hist = simulate(g0, tasks_chunk[offset], n_subs, mu, N_pop, gamma,
                        task_weights, m=m,
                        seed=make_sim_seed(rep, T, dT, m, density, alpha),
                        n_genome_snapshots=n_snaps,
                        task_sampling=task_sampling, fitness_r=fitness_r,
                        max_redraws=max_redraws,
                        record_modularity=record_modularity)
        hist['rep_index'] = rep
        results.append(hist)

    return (T, dT, m, density, alpha), rep_start, results


# ============================================================
# 13. RUNNER
# ============================================================

def ensure_task_ensembles(cache_dir, L, T, task_divs, n_ensembles):
    """Load cached ensembles, rebuilding if absent, stale, or too few."""
    tpath = task_cache_path(cache_dir, L, T)

    if os.path.exists(tpath + '.npz') and os.path.exists(tpath + '_meta.json'):
        try:
            ensembles, realized, meta = load_task_ensembles(tpath)
            have_all = all(dT in ensembles for dT in task_divs)
            if have_all and meta['n_ensembles'] >= n_ensembles:
                alpha_map = {float(k): v for k, v in meta['alpha_map'].items()}
                print(f'  Loaded task ensembles T={T} '
                      f'({meta["n_ensembles"]} per dT)')
                return alpha_map, ensembles, realized
        except (OSError, KeyError, ValueError):
            pass

    print(f'  Building {n_ensembles} task ensembles per dT for T={T}...')
    alpha_map, ensembles, realized = build_task_ensembles(
        T, L, task_divs, n_ensembles)
    save_task_ensembles(tpath, ensembles, alpha_map, realized, T, L, n_ensembles)
    for dT in task_divs:
        print(f'    dT={dT}: alpha={alpha_map[dT]:.4f}  realized='
              f'{np.mean(realized[dT]):.4f} +/- {np.std(realized[dT], ddof=1):.4f}')
    return alpha_map, ensembles, realized


def resolve_m_values(m_spec, T: int) -> List[int]:
    """Simultaneity levels for one T.

        None         every level, 1 through T
        'endpoints'  sequential and simultaneous only, {1, T}
        list         explicit levels, dropping any that exceed T
    """
    if m_spec is None:
        return list(range(1, T + 1))
    if isinstance(m_spec, str):
        if m_spec == 'endpoints':
            return sorted({1, T})
        raise ValueError(f'Unknown M_VALUES spec: {m_spec!r}')
    return [int(m) for m in m_spec if int(m) <= T]


def build_work(cfg: dict, alpha_maps: Dict, ensembles: Dict) -> List[tuple]:
    """Expand the parameter grid into chunked work items, skipping conditions
        already cached at a sufficient substitution budget.
    """
    L, K = cfg['L'], cfg['K']
    gamma, fitness_r = cfg['GAMMA'], cfg['FITNESS_R']
    n_reps, chunk = cfg['N_REPS'], cfg['CHUNK_SIZE']
    cache_dir = cfg['CACHE_DIR']
    densities = cfg['GENOME_DENSITIES'] or [1.0 / K]

    rule: Callable[[int, int], int] = cfg.get('N_SUBS_RULE') or (
        lambda T, m: (cfg['N_SUBS'] if cfg['N_SUBS'] is not None
                      else n_subs_rule(T, m, L, cfg['N_SUBS_EPOCHS'],
                                       cfg['N_SUBS_FLOOR'])))

    work = []
    for T in cfg['T_VALUES']:
        m_values = resolve_m_values(cfg['M_VALUES'], T)
        task_weights = np.ones(T) / T

        for dT in cfg['TASK_DIVERGENCES']:
            alpha = alpha_maps[T][dT]
            stack = ensembles[T][dT]

            for m in m_values:
                if m > T:
                    continue
                n_subs = int(rule(T, m))

                for density in densities:
                    spath = sim_cache_path(cache_dir, L, K, gamma, fitness_r,
                                           density, T, dT, m, alpha)
                    have = cached_condition_spec(spath)
                    if have is not None:
                        if have[0] >= n_subs and have[1] >= n_reps:
                            continue
                        print(f'Re-running (cached n_subs={have[0]}, '
                              f'n_reps={have[1]}; need {n_subs}, {n_reps}): '
                              f'T={T} dT={dT} m={m} density={density:.4f}')

                    for start in range(0, n_reps, chunk):
                        stop = min(start + chunk, n_reps)
                        work.append((
                            T, dT, m, density, alpha, stack[start:stop], start,
                            n_subs, cfg['MU'], cfg['N_POP'], gamma, K, L,
                            task_weights, cfg['N_GENOME_SNAPSHOTS'],
                            cfg['TASK_SAMPLING'], fitness_r, cfg['MAX_REDRAWS'],
                            cfg['RECORD_MODULARITY'],
                        ))
    return work


def estimate_cost(work: List[tuple]) -> Tuple[float, float]:
    """(total NNLS solves, longest single work item). A lower bound: failed
        epochs also force full enumeration and are not predictable.
    """
    per_item = []
    for w in work:
        (T, _, m, _, _, chunk_tasks, _, n_subs, _, _, _, K, L, _, n_snaps,
         _, _, _, rec_mod) = w
        n_chunk, states = chunk_tasks.shape[0], n_subs + 1
        full_states = states if rec_mod == 'all' else min(n_snaps + 2, states)
        per_item.append(n_chunk * (states * (T + L * K * m)
                                   + full_states * L * K * max(T - m, 0)))
    return float(sum(per_item)), float(max(per_item)) if per_item else 0.0


def run_simulations(cfg: dict, dry_run: bool = False):
    L, K = cfg['L'], cfg['K']
    gamma, fitness_r = cfg['GAMMA'], cfg['FITNESS_R']
    cache_dir = cfg['CACHE_DIR']
    n_workers = cfg['N_WORKERS'] or os.cpu_count()
    cfg.setdefault('RECORD_MODULARITY', DEFAULTS['RECORD_MODULARITY'])
    os.makedirs(cache_dir, exist_ok=True)

    alpha_maps, ensembles = {}, {}
    for T in cfg['T_VALUES']:
        a, e, _ = ensure_task_ensembles(cache_dir, L, T,
                                        cfg['TASK_DIVERGENCES'], cfg['N_REPS'])
        alpha_maps[T], ensembles[T] = a, e

    work = build_work(cfg, alpha_maps, ensembles)
    if not work:
        print('All conditions already cached at sufficient depth.')
        return

    total, longest = estimate_cost(work)
    n_conditions = len({(w[0], w[1], w[2], w[3]) for w in work})
    print(f'\n{len(work)} work items over {n_conditions} conditions, '
          f'{n_workers} workers')
    print(f'  total   {total / 1e9:.2f}e9 NNLS solves')
    print(f'  longest {longest / 1e6:.1f}e6 NNLS solves (critical path)')

    if dry_run:
        print('\nDry run; nothing executed.')
        for spec in sorted({(w[0], w[2], w[7]) for w in work}):
            print(f'  T={spec[0]:>2d} m={spec[1]:>2d} n_subs={spec[2]}')
        return

    expected = defaultdict(int)
    for w in work:
        expected[(w[0], w[1], w[2], w[3])] += 1
    pending: Dict[tuple, List] = defaultdict(list)
    failed = set()
    t0 = time.time()

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(_worker, a): a for a in work}

        for fut in as_completed(futures):
            spec = futures[fut]
            key4 = (spec[0], spec[1], spec[2], spec[3])
            try:
                _, rep_start, results = fut.result()
            except Exception as exc:
                print(f'  FAILED T={key4[0]} dT={key4[1]} m={key4[2]} '
                      f'density={key4[3]:.4f} rep_start={spec[6]}: {exc}')
                failed.add(key4)
                pending.pop(key4, None)   # a failed condition never saves;
                continue                  # release its buffered chunks

            if key4 in failed:
                continue
            pending[key4].append((rep_start, results))
            if len(pending[key4]) < expected[key4]:
                continue

            T, dT, m, density = key4
            alpha = spec[4]
            merged = [h for _, chunk in sorted(pending.pop(key4))
                      for h in chunk]
            merged.sort(key=lambda h: h['rep_index'])

            spath = sim_cache_path(cache_dir, L, K, gamma, fitness_r,
                                   density, T, dT, m, alpha)
            params = dict(T=T, dT=dT, m=m, density=density, alpha=alpha,
                          L=L, K=K, GAMMA=gamma, FITNESS_R=fitness_r,
                          N_SUBS=spec[7], N_REPS=cfg['N_REPS'],
                          MU=cfg['MU'], N_POP=cfg['N_POP'],
                          TASK_SEED=cfg['TASK_SEED'],
                          N_GENOME_SNAPSHOTS=cfg['N_GENOME_SNAPSHOTS'],
                          TASK_SAMPLING=cfg['TASK_SAMPLING'],
                          MAX_REDRAWS=cfg['MAX_REDRAWS'],
                          RECORD_MODULARITY=cfg['RECORD_MODULARITY'])
            save_condition(spath, merged, params)

            reasons = defaultdict(int)
            for h in merged:
                reasons[h['termination_reason']] += 1
            print(f'  Saved T={T} dT={dT} m={m} density={density:.4f}  '
                  f'{dict(reasons)}  ({time.time() - t0:.1f}s)')

    if failed:
        print(f'\n{len(failed)} conditions failed and were not saved:')
        for key4 in sorted(failed):
            print(f'  T={key4[0]} dT={key4[1]} m={key4[2]} '
                  f'density={key4[3]:.4f}')

    print(f'\nDone. Total: {time.time() - t0:.1f}s')


# ============================================================
# 14. CLI
# ============================================================

def parse_args():
    d = DEFAULTS
    p = argparse.ArgumentParser()
    p.add_argument('--L', type=int, default=d['L'])
    p.add_argument('--K', type=int, default=d['K'])
    p.add_argument('--gamma', type=float, default=d['GAMMA'])
    p.add_argument('--fitness_r', type=float, default=d['FITNESS_R'],
                   help='Power-mean exponent: 0 = geometric mean, '
                        '>0 soft-max, <0 soft-min.')
    p.add_argument('--T', type=int, nargs='+', default=d['T_VALUES'],
                   dest='T_VALUES')
    p.add_argument('--m', type=int, nargs='+', default=d['M_VALUES'],
                   dest='M_VALUES')
    p.add_argument('--m_endpoints', action='store_true',
                   help='Run only m=1 and m=T, overriding --m. Intended for '
                        'robustness specs that contrast sequential against '
                        'simultaneous selection.')
    p.add_argument('--dT', type=float, nargs='+',
                   default=d['TASK_DIVERGENCES'], dest='TASK_DIVERGENCES')
    p.add_argument('--densities', type=float, nargs='+',
                   default=d['GENOME_DENSITIES'], dest='GENOME_DENSITIES')
    p.add_argument('--n_reps', type=int, default=d['N_REPS'])
    p.add_argument('--n_pop', type=int, default=d['N_POP'])
    p.add_argument('--mu', type=float, default=d['MU'])
    p.add_argument('--n_subs', type=int, default=d['N_SUBS'],
                   help='Fixed budget for every condition. Default None uses '
                        'the T- and m-dependent rule.')
    p.add_argument('--n_subs_epochs', type=int, default=d['N_SUBS_EPOCHS'],
                   help='Target selective epochs per task; default L.')
    p.add_argument('--n_subs_floor', type=int, default=d['N_SUBS_FLOOR'])
    p.add_argument('--max_redraws', type=int, default=d['MAX_REDRAWS'])
    p.add_argument('--chunk_size', type=int, default=d['CHUNK_SIZE'])
    p.add_argument('--record_modularity', type=str,
                   default=d['RECORD_MODULARITY'],
                   choices=['snapshots', 'all', 'none'],
                   help="'all' forces full-repertoire enumeration at every "
                        'step and is substantially more expensive.')
    p.add_argument('--task_seed', type=int, default=d['TASK_SEED'])
    p.add_argument('--cache_dir', type=str, default=d['CACHE_DIR'])
    p.add_argument('--workers', type=int, default=d['N_WORKERS'])
    p.add_argument('--n_snaps', type=int, default=d['N_GENOME_SNAPSHOTS'])
    p.add_argument('--task_sampling', type=str, default=d['TASK_SAMPLING'],
                   choices=['random', 'sliding'])
    p.add_argument('--dry_run', action='store_true',
                   help='Print the work plan and cost estimate, run nothing.')
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    cfg = {
        'L': args.L, 'K': args.K, 'GAMMA': args.gamma,
        'FITNESS_R': args.fitness_r,
        'T_VALUES': args.T_VALUES,
        'M_VALUES': 'endpoints' if args.m_endpoints else args.M_VALUES,
        'TASK_DIVERGENCES': args.TASK_DIVERGENCES,
        'GENOME_DENSITIES': args.GENOME_DENSITIES,
        'N_REPS': args.n_reps, 'N_POP': args.n_pop, 'MU': args.mu,
        'N_SUBS': args.n_subs, 'N_SUBS_EPOCHS': args.n_subs_epochs,
        'N_SUBS_FLOOR': args.n_subs_floor, 'N_SUBS_RULE': None,
        'MAX_REDRAWS': args.max_redraws, 'CHUNK_SIZE': args.chunk_size,
        'RECORD_MODULARITY': args.record_modularity,
        'TASK_SEED': args.task_seed, 'CACHE_DIR': args.cache_dir,
        'N_WORKERS': args.workers, 'N_GENOME_SNAPSHOTS': args.n_snaps,
        'TASK_SAMPLING': args.task_sampling,
    }
    run_simulations(cfg, dry_run=args.dry_run)