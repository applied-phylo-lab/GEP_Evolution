#!/usr/bin/env python3
"""
figlib.py
=========
Shared analysis layer for the figure scripts.

Every figure resolves a cutoff, extracts a metric per replicate, and aggregates.
Those three steps live here so that F2, F3, F4 and the supplementary figures
cannot drift apart in how they define "after 200 substitutions" or how they
normalize differentiation. Only depends on simulate.py.

--------------------------------------------------------------------
Indexing
--------------------------------------------------------------------
State arrays are indexed by realized substitutions: index k is the genotype
after exactly k substitutions, index 0 the initial genome. A cutoff of C
therefore reads index C, not C-1.

A replicate that terminated before C is read at its last state. Under the
current engine that is exact rather than approximate: 'absorbing' means no
mutation is beneficial under any admissible active subset, so the genotype
cannot change again and the last state IS the state at every later cutoff.
The one exception is 'redraw_cap', where the walk was cut short by the
numerical safeguard rather than by reaching an optimum; those replicates are
counted separately by `termination_summary` and should be reported, not
silently pooled.

--------------------------------------------------------------------
Two cutoff currencies
--------------------------------------------------------------------
`Cutoff('substitutions', C)` compares conditions after the same number of
realized substitutions. This holds constant the amount of evolutionary change
in the shared genotype and is the primary currency for the main figures.

`Cutoff('exposure', E)` compares conditions after each task has contributed to
fitness in E selective epochs. Because a substitution scores m of T tasks,
matching E gives low-m conditions proportionally more substitutions. This is
the supplementary control for the objection that high-T, low-m populations
simply received less direct selection per task.

Exposure has two modes. 'realized' walks the cached `ep_counts` and returns the
first state at which the mean epoch count per task reaches E, so the matching is
a property of the trajectory rather than of its expectation. 'expected' uses the
closed form S = round(E*T/m) and is provided for comparison; the two agree
closely because ep_counts is multinomial around its mean, but only 'realized'
is exact.

--------------------------------------------------------------------
Normalization
--------------------------------------------------------------------
Each replicate has its own task ensemble, so differentiation is normalized by
that replicate's own realized mean pairwise task divergence
(`task_dT_realized`), not by a condition-level constant. The nominal dT is a
label for the calibration target; realized values scatter around it, and the
scatter is widest at small T where few pairs contribute to the mean.

--------------------------------------------------------------------
Pairing
--------------------------------------------------------------------
Replicate i uses the same initial genome and the same task ensemble at every m,
so contrasts across m are paired. `paired_difference` matches on `rep_index`
and returns the mean and standard error of the within-replicate differences,
which is both tighter and free of the independence assumption that an unpaired
difference-of-means would require.
"""

from dataclasses import dataclass, field
from itertools import combinations
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib as mpl
import numpy as np

from simulate import (load_condition, load_task_ensembles,
                      sim_cache_path, task_cache_path)


# ============================================================
# 1. STYLE
# ============================================================

def apply_style():
    mpl.rcParams.update({
        'pdf.use14corefonts': True,
        'font.family': 'sans-serif',
        'font.sans-serif': ['Helvetica'],
        'axes.spines.top': False,
        'axes.spines.right': False,
        'font.size': 11,
    })


def dt_colors(task_divs: Sequence[float]) -> Dict[float, tuple]:
    """One colour per task divergence, dark at high divergence."""
    vals = sorted(task_divs)
    cmap = mpl.colormaps['viridis_r']
    if len(vals) == 1:
        return {vals[0]: cmap(0.5)}
    return {dT: cmap(i / (len(vals) - 1)) for i, dT in enumerate(vals)}


def T_colors(T_values: Sequence[int]) -> Dict[int, tuple]:
    """One colour per task count. Stops short of the pale end of plasma so
    every line stays legible on white."""
    vals = sorted(T_values)
    cmap = mpl.colormaps['plasma']
    if len(vals) == 1:
        return {vals[0]: cmap(0.4)}
    return {T: cmap(0.08 + 0.72 * i / (len(vals) - 1))
            for i, T in enumerate(vals)}


def add_dt_colorbar(fig, task_divs: Sequence[float], orientation='horizontal',
                    rect=None):
    norm = mpl.colors.Normalize(vmin=min(task_divs), vmax=max(task_divs))
    cmap = mpl.colormaps['viridis_r']
    if rect is None:
        rect = ([0.27, 0.02, 0.50, 0.020] if orientation == 'horizontal'
                else [0.97, 0.15, 0.015, 0.70])
    ax = fig.add_axes(rect)
    cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm,
                                   orientation=orientation, ticks=task_divs)
    cb.set_ticklabels([f'{dT:g}' for dT in task_divs])
    if orientation == 'horizontal':
        cb.set_label(r'$\overline{\Delta T}$', fontsize=10)
    else:
        cb.ax.set_title(r'$\overline{\Delta T}$', fontsize=10, pad=4)
    cb.ax.tick_params(labelsize=8)
    return cb


# ============================================================
# 2. CACHE ACCESS
# ============================================================

@dataclass
class CacheSpec:
    """Identifies one parameter root plus a genome density."""
    cache_dir: str = 'simulation_cache'
    L: int = 100
    K: int = 4
    gamma: float = 1.0
    fitness_r: float = 0.0
    density: float = 0.25
    T_values: List[int] = field(default_factory=lambda: [2, 4, 6, 8])
    task_divs: List[float] = field(
        default_factory=lambda: [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])

    def label(self) -> str:
        return (f'L{self.L}_K{self.K}_gamma{self.gamma}'
                f'_fr{self.fitness_r}_density{self.density:.4f}')


def load_alpha_maps(spec: CacheSpec) -> Dict[int, Dict[float, float]]:
    """Calibrated Dirichlet concentration per (T, dT), read from the task
    ensemble metadata. Needed to reconstruct simulation cache filenames."""
    out = {}
    for T in spec.T_values:
        try:
            _, _, meta = load_task_ensembles(
                task_cache_path(spec.cache_dir, spec.L, T))
        except (OSError, KeyError, ValueError):
            print(f'WARNING: no task ensembles for T={T}')
            continue
        out[T] = {float(k): v for k, v in meta['alpha_map'].items()}
    return out


def load_grid(spec: CacheSpec, m_values: Optional[Sequence[int]] = None,
              verbose: bool = True) -> Dict[int, Dict[float, Dict[int, List]]]:
    """{T: {dT: {m: replicates}}} for everything present in the cache.

    Missing conditions are skipped with a warning rather than raising, so a
    partially complete batch still plots what it has.
    """
    alpha_maps = load_alpha_maps(spec)
    data: Dict[int, Dict[float, Dict[int, List]]] = {}

    for T in spec.T_values:
        if T not in alpha_maps:
            continue
        data[T] = {}
        want_m = list(range(1, T + 1)) if m_values is None else \
            [m for m in m_values if m <= T]

        for dT in spec.task_divs:
            if dT not in alpha_maps[T]:
                continue
            alpha = alpha_maps[T][dT]
            data[T][dT] = {}
            for m in want_m:
                path = sim_cache_path(spec.cache_dir, spec.L, spec.K,
                                      spec.gamma, spec.fitness_r, spec.density,
                                      T, dT, m, alpha)
                try:
                    reps, _ = load_condition(path)
                except (OSError, ValueError):
                    continue
                data[T][dT][m] = reps
            if verbose and data[T][dT]:
                print(f'  T={T} dT={dT}: m={sorted(data[T][dT])} '
                      f'({len(next(iter(data[T][dT].values())))} reps)')
    return data


def get(data: Dict, T: int, dT: float, m: int) -> Optional[List]:
    """Replicates for one condition, or None if absent."""
    return data.get(T, {}).get(dT, {}).get(m)


def resolve_m(m_selector, available: Sequence[int], T: int) -> Optional[int]:
    """'min' -> sequential, 'T' -> simultaneous, 'max' -> largest present,
    or an explicit integer."""
    available = sorted(available)
    if not available:
        return None
    if isinstance(m_selector, int):
        return m_selector if m_selector in available else None
    if m_selector == 'min':
        return min(available)
    if m_selector == 'max':
        return max(available)
    if m_selector == 'T':
        return T if T in available else None
    raise ValueError(f'Unknown m selector: {m_selector!r}')


# ============================================================
# 3. CUTOFFS
# ============================================================

@dataclass(frozen=True)
class Cutoff:
    """Where along a trajectory to read a metric. See the module docstring."""
    kind: str                      # 'substitutions' | 'exposure'
    value: int
    exposure_mode: str = 'realized'    # 'realized' | 'expected'

    def index(self, rep: Dict) -> int:
        """State index for one replicate, clamped to its last state."""
        last = int(rep['n_states']) - 1

        if self.kind == 'substitutions':
            return min(int(self.value), last)

        if self.kind == 'exposure':
            if self.exposure_mode == 'expected':
                T, m = int(rep['n_tasks']), _rep_m(rep)
                return min(int(round(self.value * T / m)), last)
            ep = np.asarray(rep['ep_counts'], dtype=float)
            reached = np.flatnonzero(ep.mean(axis=1) >= self.value)
            return int(reached[0]) if reached.size else last

        raise ValueError(f'Unknown cutoff kind: {self.kind!r}')

    def reached(self, rep: Dict) -> bool:
        """Did this replicate actually reach the requested cutoff, or was it
        clamped? Clamping is exact for 'absorbing' replicates but should be
        reported for 'redraw_cap' ones."""
        last = int(rep['n_states']) - 1
        if self.kind == 'substitutions':
            return last >= int(self.value)
        ep = np.asarray(rep['ep_counts'], dtype=float)
        return bool(ep[-1].mean() >= self.value)

    def label(self) -> str:
        if self.kind == 'substitutions':
            return f'{self.value} substitutions'
        return f'{self.value} selective epochs per task'


def _rep_m(rep: Dict) -> int:
    """Simultaneity, recovered from the width of the active-task record."""
    act = np.asarray(rep['active_tasks'])
    return int(act.shape[1]) if act.ndim == 2 and act.shape[0] else 1


# ============================================================
# 4. METRICS
# ============================================================

def differentiation(rep: Dict, cutoff: Cutoff) -> float:
    """Mean pairwise phenotype distance, normalized by this replicate's own
    realized task divergence. Bounded in [0, 1]; 1 means the evolved phenotypes
    are as far apart as the task optima they serve."""
    td = float(rep['task_dT_realized'])
    if not (np.isfinite(td) and td > 0):
        return np.nan
    return float(rep['pheno_dist'][cutoff.index(rep)]) / td


def optimization(rep: Dict, cutoff: Cutoff) -> float:
    """1 - ||d||_2 / sqrt(T), the fraction of the task deficit eliminated.
    Bounded in [0, 1] because each residual is at most 1."""
    d = np.asarray(rep['d'], dtype=float)[cutoff.index(rep)]
    return 1.0 - float(np.linalg.norm(d)) / np.sqrt(d.shape[0])


METRICS = {
    'differentiation': (differentiation, r'Degree of differentiation'),
    'optimization': (optimization, r'Degree of optimization'),
}


def metric_values(reps: Sequence[Dict], metric: str,
                  cutoff: Cutoff) -> Tuple[np.ndarray, np.ndarray]:
    """(rep_index, value) arrays for one condition. Indices are returned so
    that contrasts across conditions can be paired on replicate identity
    rather than on array position."""
    fn = METRICS[metric][0]
    idx, vals = [], []
    for rep in reps:
        v = fn(rep, cutoff)
        if np.isfinite(v):
            idx.append(int(rep['rep_index']))
            vals.append(v)
    return np.array(idx, dtype=int), np.array(vals, dtype=float)


def metric_label(metric: str) -> str:
    return METRICS[metric][1]


# ============================================================
# 5. AGGREGATION
# ============================================================

def mean_sd(vals: np.ndarray) -> Tuple[float, float]:
    """Mean and sample SD. SD describes variation in the quantity across task
    worlds, which is what the main figures plot."""
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return np.nan, np.nan
    if vals.size == 1:
        return float(vals[0]), 0.0
    return float(vals.mean()), float(vals.std(ddof=1))


def mean_se(vals: np.ndarray) -> Tuple[float, float]:
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return np.nan, np.nan
    if vals.size == 1:
        return float(vals[0]), 0.0
    return float(vals.mean()), float(vals.std(ddof=1) / np.sqrt(vals.size))


def paired_difference(idx_a: np.ndarray, vals_a: np.ndarray,
                      idx_b: np.ndarray, vals_b: np.ndarray
                      ) -> Tuple[float, float, int]:
    """(mean difference, standard error, n) for a - b, matched on rep_index.

    Replicate i shares its initial genome and task ensemble across all m, so
    the difference is paired. Pairing removes both shared components from the
    error and needs no independence assumption.
    """
    lookup = dict(zip(idx_b.tolist(), vals_b.tolist()))
    diffs = [va - lookup[i] for i, va in zip(idx_a.tolist(), vals_a.tolist())
             if i in lookup]
    diffs = np.array([d for d in diffs if np.isfinite(d)], dtype=float)
    if diffs.size == 0:
        return np.nan, np.nan, 0
    if diffs.size == 1:
        return float(diffs[0]), 0.0, 1
    return (float(diffs.mean()),
            float(diffs.std(ddof=1) / np.sqrt(diffs.size)),
            int(diffs.size))


# ============================================================
# 6. DIAGNOSTICS
# ============================================================

def termination_summary(reps: Sequence[Dict],
                        cutoff: Optional[Cutoff] = None) -> Dict:
    """Termination reasons, realized substitutions, and (if a cutoff is given)
    how many replicates were clamped rather than reaching it.

    `n_clamped` is not a defect on its own: an 'absorbing' replicate read past
    its last state is exact. It matters when `n_redraw_cap` is nonzero.
    """
    reasons: Dict[str, int] = {}
    R = []
    for rep in reps:
        reasons[rep['termination_reason']] = \
            reasons.get(rep['termination_reason'], 0) + 1
        R.append(int(rep['n_subs_realized']))
    R = np.array(R, dtype=float)

    out = {
        'n_reps': len(reps),
        'reasons': reasons,
        'n_redraw_cap': reasons.get('redraw_cap', 0),
        'R_min': float(R.min()) if R.size else np.nan,
        'R_median': float(np.median(R)) if R.size else np.nan,
        'R_max': float(R.max()) if R.size else np.nan,
    }
    if cutoff is not None:
        out['n_clamped'] = int(sum(1 for r in reps if not cutoff.reached(r)))
    return out


def exposure_summary(reps: Sequence[Dict], cutoff: Cutoff) -> Dict:
    """Realized selective epochs per task at the cutoff: the mean across tasks
    and the minimum over tasks, both averaged over replicates. Under random
    task sampling the per-task counts are multinomial around the mean, so the
    minimum is the honest check on whether matching held."""
    means, mins = [], []
    for rep in reps:
        ep = np.asarray(rep['ep_counts'], dtype=float)[cutoff.index(rep)]
        means.append(ep.mean())
        mins.append(ep.min())
    if not means:
        return {'ep_mean': np.nan, 'ep_min': np.nan}
    return {'ep_mean': float(np.mean(means)), 'ep_min': float(np.mean(mins))}


def realized_dT_summary(reps: Sequence[Dict]) -> Tuple[float, float]:
    """Mean and SD of the realized task divergence across replicates. The
    nominal dT is a calibration target; this is what the ensembles actually
    delivered."""
    return mean_sd(np.array([rep['task_dT_realized'] for rep in reps]))


# ============================================================
# 7. MISC
# ============================================================

def mean_pairwise_distance(X: np.ndarray) -> float:
    """Mean pairwise Euclidean distance between the columns of X. Retained for
    ad-hoc checks; replicates carry their own `task_dT_realized`."""
    T = X.shape[1]
    if T < 2:
        return np.nan
    return float(np.mean([np.linalg.norm(X[:, i] - X[:, j])
                          for i, j in combinations(range(T), 2)]))


def panel_label(i: int) -> str:
    """A, B, ... Z, AA, AB, ... so a figure never runs out of panel labels."""
    s = ''
    i += 1
    while i > 0:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s
