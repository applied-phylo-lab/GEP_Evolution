#!/usr/bin/env python3
"""
simulate.py
===========
SSWM simulation of gene expression program evolution under variable
simultaneity of selection.

Global parameters (all tunable via CLI):
  L                  : number of molecular traits
  K                  : number of gene expression programs
  GAMMA              : performance exponent P = (1-d)^gamma
  T_VALUES           : task counts to sweep
  M_VALUES           : simultaneity levels to sweep (1 <= m <= T);
                       None -> [1, 2, ..., T] for each T
  TASK_DIVERGENCES   : target mean pairwise distances between tasks
  TASK_ALPHAS        : explicit Dirichlet alphas (same length as TASK_DIVERGENCES);
                       None -> calibrate alpha to hit each target dT
  GENOME_DENSITIES   : initial genome densities to sweep; None -> [1/K]
  N_SUBS             : substitution steps per replicate
  N_REPS             : replicates per condition
  N_POP              : effective population size
  MU                 : mutation rate
  TASK_SEED          : RNG seed for task generation (fixed across runs)
  N_GENOME_SNAPSHOTS : number of intermediate snapshots (excl. initial/final)

Cache layout:
  simulation_cache/
    L{L}_K{K}_gamma{GAMMA}_fr{FITNESS_R}/
      tasks_T{T}.npz                             <- task vectors, once per T
      tasks_T{T}_meta.json
      density{DENSITY}/
        sim_T{T}_dT{dT}_m{m}_alpha{ALPHA}.npz   <- all reps for one condition
        sim_T{T}_dT{dT}_m{m}_alpha{ALPHA}_meta.json

Snapshot structure (identical for initial, intermediate, final timepoints):
  genome   : (L, K) binary array
  usage    : (T, K) activation vectors, one row per task
  phenotype: (L, T) expressed phenotypes G @ a_j, one column per task
  cum_time : float, cumulative evolutionary time at the moment this genome
             was the current state (i.e. time of its fixation)

Usage:
  python simulate.py
  python simulate.py --L 100 --K 6 --T 3 6 9 --m 1 2 3 6 --densities 0.167 0.333
  python simulate.py --alphas 0.01 0.01 0.01   # sparse task robustness check
"""

import argparse
import hashlib
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.optimize import nnls


# ============================================================
# 0. SEED GENERATION
# ============================================================

DEFAULTS = dict(
    N_POP=10_000,
    MU=1e-7,
    L=100,
    K=4,
    GAMMA=1.0,
    N_SUBS=500,
    N_REPS=50,
    T_VALUES=[2, 3, 4, 6, 8],
    M_VALUES=None,           # None -> [1, 2, ..., T] for each T
    TASK_DIVERGENCES=[0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4],
    TASK_ALPHAS=None,        # None -> calibrate alpha to hit TASK_DIVERGENCES
                             # explicit list -> use directly, bypasses dT targeting
                             # useful for sparse task robustness checks
    GENOME_DENSITIES=[0.25],   # None -> [1/K]
    TASK_SEED=270,
    CACHE_DIR='simulation_cache',
    N_WORKERS=100,          # None -> os.cpu_count()
    N_GENOME_SNAPSHOTS=10,   # intermediate snapshots (excl. initial and final)
    TASK_SAMPLING='random',  # 'random' or 'sliding' (deterministic cyclic window)
    FITNESS_R=0.0,           # power mean exponent: 0.0 = geometric mean (limit),
                             # >0 = soft-max direction, <0 = soft-min direction
)


# ============================================================
# 2. TASK GENERATION
# ============================================================

def create_dirichlet_tasks(L: int, T: int, alpha: float, seed: int) -> np.ndarray:
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


def build_task_map(T: int, L: int, task_divergences: List[float],
                   task_seed: int,
                   task_alphas: Optional[List[float]] = None,
                   ) -> Tuple[Dict, Dict]:
    """
    Build task map for a given T.

    If task_alphas is None, alpha is calibrated to hit each target dT.
    If task_alphas is provided explicitly, it must have the same length as
    task_divergences — the dT values are then used only as labels/keys,
    and the actual alpha values are used directly. This allows sparse task
    robustness checks by setting low alpha independently of dT targeting.
    """
    if task_alphas is not None:
        assert len(task_alphas) == len(task_divergences), \
            "task_alphas must have same length as task_divergences"
        alpha_map = dict(zip(task_divergences, task_alphas))
    else:
        alpha_map = {dT: find_alpha_for_target(dT, L, T)
                     for dT in task_divergences}

    task_map = {dT: create_dirichlet_tasks(L, T, alpha_map[dT], seed=task_seed)
                for dT in task_divergences}
    return alpha_map, task_map


# ============================================================
# 3. GENOME
# ============================================================

def make_genome(L: int, K: int, density: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return (rng.random((L, K)) < density).astype(float)

def make_genome_seed(L: int, K: int, density: float) -> int:
    key = f'genome_{L}_{K}_{density:.10g}'
    return int(hashlib.md5(key.encode()).hexdigest(), 16) % (2 ** 31)


def make_sim_seed(rep: int, T: int, dT: float, m: int,
                  density: float, alpha: float) -> int:
    key = f'sim_{rep}_{T}_{dT:.10g}_{m}_{density:.10g}_{alpha:.10g}'
    return int(hashlib.md5(key.encode()).hexdigest(), 16) % (2 ** 31)

# ============================================================
# 4. PERFORMANCE
# ============================================================

def compute_performance(genome: np.ndarray, tasks: np.ndarray,
                        gamma: float) -> Dict:
    """
    Returns dict with:
        d        : (T,) residual distances
        P        : (T,) performances
        a        : list of T activation vectors, each (K,)
        usage    : (T, K) stacked activation vectors
        phenotype: (L, T) expressed phenotypes G @ a_j
    """
    n_tasks = tasks.shape[1]
    d = np.zeros(n_tasks)
    P = np.zeros(n_tasks)
    a_list = []
    phenotypes = np.zeros((genome.shape[0], n_tasks))

    for e in range(n_tasks):
        a, _ = nnls(genome, tasks[:, e])
        z = genome @ a
        d[e] = np.linalg.norm(tasks[:, e] - z)
        P[e] = max(0.0, 1.0 - d[e]) ** gamma
        a_list.append(a)
        phenotypes[:, e] = z

    return {
        'd':         d,
        'P':         P,
        'a':         a_list,
        'usage':     np.array(a_list),  # (T, K)
        'phenotype': phenotypes,        # (L, T)
    }


def make_snapshot(genome: np.ndarray, tasks: np.ndarray, gamma: float) -> Dict:
    """
    Identical snapshot structure for initial, intermediate, and final timepoints.
    """
    info = compute_performance(genome, tasks, gamma)
    return {
        'genome':    genome.copy(),
        'usage':     info['usage'],
        'phenotype': info['phenotype'],
    }


# ============================================================
# 5. FITNESS
# ============================================================

def fitness_power_mean(P: np.ndarray, w: np.ndarray, r: float) -> float:
    """
    Weighted power mean of task performances with exponent r.

    r = 0.0  : geometric mean (exact limit of power mean as r -> 0).
               Returns 0 if any P <= 0.
    r > 0    : soft-max direction; larger r weights the best task more.
    r < 0    : soft-min direction; more negative r weights the worst task more.
               Returns 0 if any P <= 0 (log-domain protection for negative r).

    For any r, returns 0 if P contains non-positive values and r <= 0,
    since those cases are either undefined (log) or yield 0 (0^r for r<0 -> inf).
    """
    if r == 0.0:
        if np.any(P <= 0):
            return 0.0
        return float(np.exp(np.sum(w * np.log(P))))
    else:
        if np.any(P <= 0) and r < 0:
            return 0.0
        return float(np.power(np.sum(w * np.power(P, r)), 1.0 / r))


# ============================================================
# 6. POPULATION GENETICS
# ============================================================

def selection_coeff(W_mut: float, W_wt: float) -> float:
    """Relative selection coefficient s = (W_mut - W_wt) / W_wt."""
    if W_wt <= 0:
        return float('inf') if W_mut > 0 else 0.0
    return (W_mut - W_wt) / W_wt


def p_fix(s: float, N: int) -> float:
    """
    Fixation probability of a beneficial mutation under SSWM.
    Uses the standard Kimura formula: (1 - exp(-2s)) / (1 - exp(-2Ns)).
    Handles limiting cases for numerical stability.
    """
    if s <= 0:
        return 0.0
    Ns = N * s
    if Ns > 100:
        return 1.0 - np.exp(-2.0 * s)
    if Ns < 1e-6:
        return (1.0 / N) * (1.0 + (N - 1) * s)
    return (1.0 - np.exp(-2.0 * s)) / (1.0 - np.exp(-2.0 * Ns))


# ============================================================
# 7. METRICS
# ============================================================

def effective_rank(a_list: List[np.ndarray]) -> float:
    """
    Effective rank of the task activation matrix.

    Measures how many independent dimensions the genome uses when deploying
    its K programs across T tasks. Computed as the exponentiated Shannon
    entropy of the squared singular values of the row-normalized activation
    matrix A (T x K), where each row is the L2-normalized activation vector
    for one task.

    Returns 1.0 if all tasks recruit the same program combination,
    and min(T, K) if each task recruits a fully independent combination.
    """
    if len(a_list) < 2:
        return 1.0
    A = np.array(a_list)
    norms = np.linalg.norm(A, axis=1, keepdims=True)
    norms = np.where(norms < 1e-12, 1.0, norms)
    A_norm = A / norms
    try:
        _, s, _ = np.linalg.svd(A_norm, full_matrices=False)
    except np.linalg.LinAlgError:
        return float('nan')
    s_sq = s ** 2
    total = s_sq.sum()
    if total < 1e-12:
        return 1.0
    p = s_sq / total
    p = p[p > 1e-12]
    return float(np.exp(-np.sum(p * np.log(p))))


def modularity_entropy(dF: np.ndarray) -> float:
    """
    Mutational modularity of the current genome.

    For each possible single-bit mutation, computes how its fitness effects
    are distributed across T tasks. A mutation that affects only one task
    has entropy 0; one that affects all tasks equally has entropy log(T).
    Modularity M = 1 - mean(H) / log(T), so M=1 means all mutations are
    fully task-specific and M=0 means all mutations affect all tasks equally.

    dF : (L*K, T) array of fitness effect vectors, one row per mutation.
    Only mutations with nonzero total effect are included in the average.
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
    """
    Mean pairwise Euclidean distance between expressed phenotypes.
    phenotype : (L, T) array — one column per task.
    Measures how distinct the realized phenotypes are from each other.
    """
    T = phenotype.shape[1]
    if T < 2:
        return 0.0
    dists = [np.linalg.norm(phenotype[:, i] - phenotype[:, j])
             for i, j in combinations(range(T), 2)]
    return float(np.mean(dists))


# ============================================================
# 8. CORE SSWM ENGINE
# ============================================================

def _run_sswm(
    genome_init: np.ndarray,
    tasks: np.ndarray,
    n_subs: int,
    mu: float,
    N: int,
    gamma: float,
    m: int,
    w: np.ndarray,
    seed: int,
    n_genome_snapshots: int = 10,
    task_sampling: str = 'random',
    fitness_r: float = 0.0,
) -> Dict:
    """
    Unified SSWM engine for all simultaneity levels.
    m == T  : full simultaneity (multicellular)
    m == 1  : pure sequential (unicellular)
    1 < m < T: partial simultaneity

    At each step, m tasks are selected. Task selection mode:
      'random'  : m tasks drawn uniformly without replacement (default)
      'sliding' : deterministic cyclic window of m consecutive tasks,
                  advancing by 1 each step. Tasks are arranged in a
                  cycle [0, 1, ..., T-1, 0, 1, ...].

    Fitness = weighted power mean of task performances with exponent fitness_r.
      fitness_r = 0.0 : geometric mean (default)
      fitness_r > 0   : soft-max direction
      fitness_r < 0   : soft-min direction
    """
    rng = np.random.default_rng(seed)
    L, K = genome_init.shape
    n_tasks = tasks.shape[1]

    # Snapshot scheduling: thresholds are evenly spaced over [0, n_subs),
    # so snapshots are taken at steps 0, ~n_subs/(n+1), ~2*n_subs/(n+1), ...
    # Step 0 is always captured (threshold starts at 0.0).
    # Replicates that stall early will only hit the thresholds that fall
    # within their actual run length; the final step is always guaranteed
    # by a post-loop check regardless of when the replicate stops.
    snap_interval = n_subs / (n_genome_snapshots + 1)
    next_snap_threshold = 0.0

    hist = {
        'eff_rank':           [],
        'modularity_entropy': [],
        'pheno_dist':         [],
        'cum_time':           [],
        'wait_time':          [],
        'P':                  [],
        'd':                  [],
        'W':                  [],
        's_max':              [],
        'n_ben':              [],
        'active_tasks':       [],
    }

    snapshots: Dict[int, Dict] = {}
    cumtime = 0.0
    genome = genome_init.copy()
    actual_subs = 0

    for step in range(n_subs):
        info_wt = compute_performance(genome, tasks, gamma)
        P_wt = info_wt['P'].copy()

        # draw active task subset
        if m == n_tasks:
            active = np.arange(n_tasks)
        elif task_sampling == 'sliding':
            start = step % n_tasks
            active = np.array([i % n_tasks for i in range(start, start + m)])
        else:  # 'random'
            active = rng.choice(n_tasks, size=m, replace=False)

        w_active = w[active] / w[active].sum()
        W_wt = fitness_power_mean(P_wt[active], w_active, fitness_r)

        # record trajectory
        hist['eff_rank'].append(effective_rank(info_wt['a']))
        hist['P'].append(P_wt.copy())
        hist['d'].append(info_wt['d'].copy())
        hist['W'].append(W_wt)
        hist['cum_time'].append(cumtime)
        hist['active_tasks'].append(active.copy())
        hist['pheno_dist'].append(
            mean_pairwise_phenotype_distance(info_wt['phenotype']))

        # snapshot: trigger when step crosses next evenly-spaced threshold
        if step >= next_snap_threshold:
            snap = make_snapshot(genome, tasks, gamma)
            snap['cum_time'] = cumtime
            snapshots[step] = snap
            next_snap_threshold += snap_interval

        # evaluate all L*K mutations
        dF_all = np.zeros((L * K, n_tasks))
        beneficial = []
        idx = 0
        for i in range(L):
            for j in range(K):
                mut = genome.copy()
                mut[i, j] = 1.0 - mut[i, j]
                info_mut = compute_performance(mut, tasks, gamma)
                dF_all[idx] = info_mut['P'] - P_wt
                W_mut = fitness_power_mean(info_mut['P'][active], w_active, fitness_r)
                s = selection_coeff(W_mut, W_wt)
                if s > 0:
                    beneficial.append((i, j, s, p_fix(s, N), mut))
                idx += 1

        hist['modularity_entropy'].append(modularity_entropy(dF_all))
        hist['s_max'].append(max((b[2] for b in beneficial), default=0.0))
        hist['n_ben'].append(len(beneficial))

        # stall check
        if len(beneficial) == 0:
            hist['wait_time'].append(float('inf'))
            for _ in range(step + 1, n_subs):
                for key in hist:
                    if hist[key]:
                        hist[key].append(hist[key][-1])
            actual_subs = step + 1
            break

        # substitution
        pfix_vals = np.array([b[3] for b in beneficial])
        lambda_tot = N * mu * pfix_vals.sum()
        wait = rng.exponential(1.0 / lambda_tot)
        hist['wait_time'].append(wait)
        cumtime += wait

        probs = pfix_vals / pfix_vals.sum()
        chosen = rng.choice(len(beneficial), p=probs)
        genome = beneficial[chosen][4]
        actual_subs = step + 1

    else:
        actual_subs = n_subs

    # ensure final snapshot is always saved
    final_step = actual_subs - 1
    if final_step not in snapshots:
        snap = make_snapshot(genome, tasks, gamma)
        snap['cum_time'] = cumtime
        snapshots[final_step] = snap

    # convert lists to arrays
    for key in hist:
        hist[key] = np.array(hist[key])

    hist['snapshots']     = snapshots
    hist['n_actual_subs'] = actual_subs
    hist['n_tasks']       = n_tasks

    return hist


# ============================================================
# 9. PUBLIC WRAPPER
# ============================================================

def simulate(genome_init: np.ndarray, tasks: np.ndarray,
             n_subs: int, mu: float, N: int, gamma: float,
             task_weights: np.ndarray, m: int, seed: int,
             n_genome_snapshots: int = 10,
             task_sampling: str = 'random',
             fitness_r: float = 0.0) -> Dict:
    """
    Public entry point for all simultaneity regimes.

    m=1    : purely sequential — one task selected per substitution
    m=T    : fully simultaneous — all tasks evaluated jointly
    1<m<T  : partial simultaneity — m tasks selected per step

    task_sampling:
      'random'  : m tasks drawn uniformly without replacement (default)
      'sliding' : deterministic cyclic window of m consecutive tasks

    fitness_r: power mean exponent for fitness aggregation over tasks.
      0.0  -> geometric mean (default)
      >0   -> soft-max direction
      <0   -> soft-min direction

    task_weights are normalized internally; uniform weights are appropriate
    for the standard sweep.
    """
    w = task_weights / task_weights.sum()
    return _run_sswm(
        genome_init, tasks, n_subs, mu, N, gamma,
        m=m, w=w, seed=seed,
        n_genome_snapshots=n_genome_snapshots,
        task_sampling=task_sampling,
        fitness_r=fitness_r,
    )


# ============================================================
# 10. CACHE HELPERS
# ============================================================

def _param_root(cache_dir: str, L: int, K: int, gamma: float,
                fitness_r: float) -> str:
    folder = os.path.join(cache_dir, f'L{L}_K{K}_gamma{gamma}_fr{fitness_r}')
    os.makedirs(folder, exist_ok=True)
    return folder


def _density_root(cache_dir: str, L: int, K: int, gamma: float,
                  fitness_r: float, density: float) -> str:
    folder = os.path.join(_param_root(cache_dir, L, K, gamma, fitness_r),
                          f'density{density:.4f}')
    os.makedirs(folder, exist_ok=True)
    return folder


def task_cache_path(cache_dir: str, L: int, K: int, gamma: float,
                    fitness_r: float, T: int) -> str:
    return os.path.join(_param_root(cache_dir, L, K, gamma, fitness_r),
                        f'tasks_T{T}')


def sim_cache_path(cache_dir: str, L: int, K: int, gamma: float,
                   fitness_r: float, density: float, T: int, dT: float,
                   m: int, alpha: float) -> str:
    return os.path.join(
        _density_root(cache_dir, L, K, gamma, fitness_r, density),
        f'sim_T{T}_dT{dT}_m{m}_alpha{alpha:.4f}'
    )


def save_task_map(base_path: str, task_map: Dict, alpha_map: Dict,
                  T: int, L: int, task_seed: int):
    """Save task vectors as npz and alpha/metadata as json. base_path has no extension."""
    arrays = {f'tasks_dT{dT}': arr for dT, arr in task_map.items()}
    np.savez_compressed(base_path + '.npz', **arrays)
    meta = {
        'T': T, 'L': L, 'task_seed': task_seed,
        'dT_values': list(task_map.keys()),
        'alpha_map': {str(k): v for k, v in alpha_map.items()},
    }
    with open(base_path + '_meta.json', 'w') as f:
        json.dump(meta, f, indent=2)


def load_task_map(base_path: str) -> Dict:
    """Load task vectors from npz. Returns dict keyed by dT value."""
    data = np.load(base_path + '.npz')
    with open(base_path + '_meta.json') as f:
        meta = json.load(f)
    return {dT: data[f'tasks_dT{dT}'] for dT in meta['dT_values']}


def _flatten_snapshots(snapshots: Dict[int, Dict]) -> Dict[str, np.ndarray]:
    """
    Convert snapshots dict to flat key-value pairs for npz storage.
    Keys follow the pattern: snap_{step}_{field}
    """
    flat = {}
    for step, snap in snapshots.items():
        flat[f'snap_{step}_genome']    = snap['genome']
        flat[f'snap_{step}_usage']     = snap['usage']
        flat[f'snap_{step}_phenotype'] = snap['phenotype']
        flat[f'snap_{step}_cum_time']  = np.array([snap['cum_time']])
    return flat


def _unflatten_snapshots(arrays: Dict, step_keys: List[int]) -> Dict[int, Dict]:
    """Reconstruct snapshots dict from flat npz arrays given known step indices."""
    return {
        step: {
            'genome':    arrays[f'snap_{step}_genome'],
            'usage':     arrays[f'snap_{step}_usage'],
            'phenotype': arrays[f'snap_{step}_phenotype'],
            'cum_time':  float(arrays[f'snap_{step}_cum_time'][0]),
        }
        for step in step_keys
    }


def save_condition(base_path: str, results: List[Dict], params: Dict):
    """
    Save all replicates for one (T, dT, m, density, alpha) condition.
    Trajectory arrays and flattened snapshots go to npz.
    Metadata and per-replicate scalars go to json.
    base_path has no extension.
    """
    arrays = {}
    rep_meta = []

    for rep_idx, hist in enumerate(results):
        prefix = f'rep{rep_idx}_'

        for key in ('eff_rank', 'modularity_entropy', 'pheno_dist',
                    'cum_time', 'wait_time', 'P', 'd', 'W',
                    's_max', 'n_ben', 'active_tasks'):
            arrays[prefix + key] = hist[key]

        for k, v in _flatten_snapshots(hist['snapshots']).items():
            arrays[prefix + k] = v

        rep_meta.append({
            'n_actual_subs':  int(hist['n_actual_subs']),
            'n_tasks':        int(hist['n_tasks']),
            'snapshot_steps': sorted(hist['snapshots'].keys()),
        })

    np.savez_compressed(base_path + '.npz', **arrays)
    with open(base_path + '_meta.json', 'w') as f:
        json.dump({'params': params, 'n_reps': len(results),
                   'rep_meta': rep_meta}, f, indent=2)


def load_condition(base_path: str) -> Tuple[List[Dict], Dict]:
    """
    Load a saved condition. Returns (results, params).
    results is a list of N_REPS hist dicts with the same structure
    produced by _run_sswm, including snapshots keyed by step index.
    """
    arrays = dict(np.load(base_path + '.npz', allow_pickle=False))
    with open(base_path + '_meta.json') as f:
        meta = json.load(f)

    results = []
    for rep_idx in range(meta['n_reps']):
        prefix   = f'rep{rep_idx}_'
        rep_info = meta['rep_meta'][rep_idx]

        hist = {key: arrays[prefix + key]
                for key in ('eff_rank', 'modularity_entropy', 'pheno_dist',
                            'cum_time', 'wait_time', 'P', 'd', 'W',
                            's_max', 'n_ben', 'active_tasks')}

        snap_arrays = {k[len(prefix):]: v
                       for k, v in arrays.items()
                       if k.startswith(prefix + 'snap_')}
        hist['snapshots']     = _unflatten_snapshots(
            snap_arrays, rep_info['snapshot_steps'])
        hist['n_actual_subs'] = rep_info['n_actual_subs']
        hist['n_tasks']       = rep_info['n_tasks']
        results.append(hist)

    return results, meta['params']


# ============================================================
# 11. WORKER
# ============================================================

def _worker(args: tuple):
    """
    Process pool worker. Runs N_REPS replicates for one
    (T, dT, m, density, alpha) condition and returns results.
    Called by ProcessPoolExecutor in run_simulations.

    Initial genomes are controlled across conditions:
    for a given (density, L, K), the same starting genome is used
    regardless of T, dT, m, or alpha.

    Simulation seeds are condition- and replicate-specific.
    """
    (T, dT, m, density, alpha, tasks,
     n_subs, n_reps, mu, N_pop, gamma,
     K, L, task_weights, n_genome_snapshots,
     task_sampling, fitness_r) = args

    results = []
    for rep in range(n_reps):
        g0 = make_genome(
            L, K, density,
            seed=make_genome_seed(L, K, density)
        )

        hist = simulate(
            g0, tasks, n_subs, mu, N_pop, gamma,
            task_weights, m=m,
            seed=make_sim_seed(rep, T, dT, m, density, alpha),
            n_genome_snapshots=n_genome_snapshots,
            task_sampling=task_sampling,
            fitness_r=fitness_r,
        )
        results.append(hist)

    return T, dT, m, density, alpha, results

# ============================================================
# 12. MAIN RUNNER
# ============================================================

def run_simulations(cfg: dict):
    L         = cfg['L']
    K         = cfg['K']
    gamma     = cfg['GAMMA']
    fitness_r = cfg.get('FITNESS_R', 0.0)
    n_subs    = cfg['N_SUBS']
    n_reps    = cfg['N_REPS']
    mu        = cfg['MU']
    N_pop     = cfg['N_POP']
    task_divs = cfg['TASK_DIVERGENCES']
    task_alphas = cfg.get('TASK_ALPHAS', None)  # None -> calibrate from dT
    T_values  = cfg['T_VALUES']
    task_seed = cfg['TASK_SEED']
    cache_dir = cfg['CACHE_DIR']
    n_workers = cfg['N_WORKERS'] or os.cpu_count()
    n_snaps   = cfg['N_GENOME_SNAPSHOTS']
    densities = cfg['GENOME_DENSITIES'] or [1.0 / K]
    task_sampling = cfg.get('TASK_SAMPLING', 'random')

    os.makedirs(cache_dir, exist_ok=True)

    # build / load task maps
    task_maps  = {}
    alpha_maps = {}
    for T in T_values:
        tpath = task_cache_path(cache_dir, L, K, gamma, fitness_r, T)
        if os.path.exists(tpath + '.npz') and task_alphas is None:
            task_maps[T] = load_task_map(tpath)
            with open(tpath + '_meta.json') as f:
                meta = json.load(f)
            alpha_maps[T] = {float(k): v
                             for k, v in meta['alpha_map'].items()}
            print(f'Loaded task map T={T}')
        else:
            print(f'Building task map T={T}...')
            alpha_maps[T], task_maps[T] = build_task_map(
                T, L, task_divs, task_seed,
                task_alphas=task_alphas)
            save_task_map(tpath, task_maps[T], alpha_maps[T],
                          T, L, task_seed)
            print(f'  Saved task map T={T}')

    # build work list
    work = []

    for T in T_values:
        m_values     = cfg['M_VALUES'] if cfg['M_VALUES'] else list(range(1, T + 1))
        task_weights = np.ones(T) / T

        for dT in task_divs:
            tasks = task_maps[T][dT]
            alpha = alpha_maps[T][dT]

            for m in m_values:
                if m > T:
                    continue

                for density in densities:
                    spath = sim_cache_path(
                        cache_dir, L, K, gamma, fitness_r, density, T, dT, m, alpha)
                    if os.path.exists(spath + '.npz'):
                        print(f'Skip (cached): T={T} dT={dT} m={m} '
                              f'alpha={alpha:.4f} density={density:.4f}')
                        continue

                    work.append((
                        T, dT, m, density, alpha, tasks,
                        n_subs, n_reps, mu, N_pop, gamma,
                        K, L, task_weights, n_snaps,
                        task_sampling, fitness_r,
                    ))

    if not work:
        print('All conditions already cached.')
        return

    print(f'\nRunning {len(work)} conditions '
          f'with {n_workers} workers...')
    t0 = time.time()

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(_worker, a): a for a in work}

        for fut in as_completed(futures):
            try:
                T, dT, m, density, alpha, results = fut.result()
            except Exception as exc:
                args_failed = futures[fut]
                raise RuntimeError(
                    f"Worker crashed: T={args_failed[0]} dT={args_failed[1]} "
                    f"m={args_failed[2]} density={args_failed[3]:.4f} "
                    f"alpha={args_failed[4]:.4f}"
                ) from exc
            spath = sim_cache_path(
                cache_dir, L, K, gamma, fitness_r, density, T, dT, m, alpha)
            params = dict(T=T, dT=dT, m=m, density=density, alpha=alpha,
                          L=L, K=K, GAMMA=gamma,
                          FITNESS_R=fitness_r,
                          N_SUBS=n_subs, N_REPS=n_reps,
                          MU=mu, N_POP=N_pop,
                          TASK_SEED=task_seed,
                          N_GENOME_SNAPSHOTS=n_snaps,
                          TASK_SAMPLING=task_sampling)
            save_condition(spath, results, params)
            print(f'  Saved T={T} dT={dT} m={m} '
                  f'alpha={alpha:.4f} density={density:.4f}  '
                  f'({time.time()-t0:.1f}s)')

    print(f'\nDone. Total: {time.time()-t0:.1f}s')


# ============================================================
# 13. CLI
# ============================================================

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--L',         type=int,   default=DEFAULTS['L'])
    p.add_argument('--K',         type=int,   default=DEFAULTS['K'])
    p.add_argument('--gamma',     type=float, default=DEFAULTS['GAMMA'])
    p.add_argument('--T',         type=int,   nargs='+',
                   default=DEFAULTS['T_VALUES'],        dest='T_VALUES')
    p.add_argument('--m',         type=int,   nargs='+',
                   default=DEFAULTS['M_VALUES'],        dest='M_VALUES')
    p.add_argument('--dT',        type=float, nargs='+',
                   default=DEFAULTS['TASK_DIVERGENCES'],dest='TASK_DIVERGENCES')
    p.add_argument('--alphas',    type=float, nargs='+',
                   default=DEFAULTS['TASK_ALPHAS'],     dest='TASK_ALPHAS',
                   help='Explicit Dirichlet alphas (same length as --dT). '
                        'If provided, bypasses dT-based calibration. '
                        'Useful for sparse task robustness checks.')
    p.add_argument('--densities', type=float, nargs='+',
                   default=DEFAULTS['GENOME_DENSITIES'],dest='GENOME_DENSITIES')
    p.add_argument('--n_subs',    type=int,   default=DEFAULTS['N_SUBS'])
    p.add_argument('--n_reps',    type=int,   default=DEFAULTS['N_REPS'])
    p.add_argument('--n_pop',     type=int,   default=DEFAULTS['N_POP'])
    p.add_argument('--mu',        type=float, default=DEFAULTS['MU'])
    p.add_argument('--task_seed', type=int,   default=DEFAULTS['TASK_SEED'])
    p.add_argument('--cache_dir', type=str,   default=DEFAULTS['CACHE_DIR'])
    p.add_argument('--workers',   type=int,   default=DEFAULTS['N_WORKERS'])
    p.add_argument('--n_snaps',   type=int,
                   default=DEFAULTS['N_GENOME_SNAPSHOTS'])
    p.add_argument('--task_sampling', type=str,
                   default=DEFAULTS['TASK_SAMPLING'],
                   choices=['random', 'sliding'],
                   help="Task selection mode: 'random' (uniform draw) "
                        "or 'sliding' (deterministic cyclic window).")
    p.add_argument('--fitness_r', type=float, default=DEFAULTS['FITNESS_R'],
                   help='Power mean exponent for fitness aggregation over tasks. '
                        '0.0 = geometric mean (default); '
                        '>0 = soft-max direction; <0 = soft-min direction.')
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    cfg = {
        'L':                  args.L,
        'K':                  args.K,
        'GAMMA':              args.gamma,
        'FITNESS_R':          args.fitness_r,
        'T_VALUES':           args.T_VALUES,
        'M_VALUES':           args.M_VALUES,
        'TASK_DIVERGENCES':   args.TASK_DIVERGENCES,
        'TASK_ALPHAS':        args.TASK_ALPHAS,
        'GENOME_DENSITIES':   args.GENOME_DENSITIES,
        'N_SUBS':             args.n_subs,
        'N_REPS':             args.n_reps,
        'N_POP':              args.n_pop,
        'MU':                 args.mu,
        'TASK_SEED':          args.task_seed,
        'CACHE_DIR':          args.cache_dir,
        'N_WORKERS':          args.workers,
        'N_GENOME_SNAPSHOTS': args.n_snaps,
        'TASK_SAMPLING':      args.task_sampling,
    }
    run_simulations(cfg)