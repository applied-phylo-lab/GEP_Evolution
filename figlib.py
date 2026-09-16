#!/usr/bin/env python3
"""
figlib.py
=========
Shared analysis and plotting layer for the figure scripts. Cutoffs, metrics and
aggregation live here so that the figures cannot drift apart in how "after N
substitutions" is resolved or how differentiation is normalized.

Indexing. State index k is the genotype after k substitutions, so a cutoff of C
reads index C, not C-1. A replicate that terminated before C is read at its last
state; see `termination_summary`.

Cutoff currencies. 'substitutions' compares conditions after the same amount of
change in the shared genotype and is the primary currency. 'exposure' compares
after each task has contributed to fitness in E selective epochs.

Normalization. Differentiation is divided by each replicate's own realized task
divergence, not by the nominal dT.

Pairing. Replicate i shares its initial genome and task ensemble across every m,
so `paired_difference` matches on rep_index.

Visual grammar. Colour is task divergence, linestyle distinguishes series that
share a panel, marker is program number, and dispersion is a shaded band.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib as mpl
import numpy as np

from simulate import (load_condition, load_task_ensembles,
                      sim_cache_path, task_cache_path)


# ============================================================
# 1. STYLE
# ============================================================

# Helvetica on a Mac, its metric clone elsewhere. Set explicitly so that text
# and math resolve to the same face.
_SANS_PREFERENCE = ('Helvetica', 'Nimbus Sans', 'Arial', 'Liberation Sans',
                    'Arimo', 'DejaVu Sans')


def _is_bold(entry) -> bool:
    w = entry.weight
    return w in ('bold', 'heavy', 'black') or (isinstance(w, (int, float)) and w >= 600)


def sans_family() -> str:
    """First available Helvetica-metric family on this machine that registers a
        bold and an italic face. A family offering only its regular face is
        skipped: matplotlib silently renders bold titles and italic math at
        regular weight, which changes the figure without failing.
    """
    from matplotlib import font_manager
    faces = {}
    for f in font_manager.fontManager.ttflist:
        has_bold, has_italic = faces.get(f.name, (False, False))
        faces[f.name] = (has_bold or _is_bold(f),
                         has_italic or f.style in ('italic', 'oblique'))
    for name in _SANS_PREFERENCE:
        if all(faces.get(name, (False, False))):
            return name
    return 'DejaVu Sans'


def apply_style():
    fam = sans_family()
    mpl.rcParams.update({
        'pdf.use14corefonts': False,   # embed, as most journals require
        'font.family': 'sans-serif',
        'font.sans-serif': [fam] + [f for f in _SANS_PREFERENCE if f != fam],
        # math must resolve to the same face as the surrounding text
        'mathtext.fontset': 'custom',
        'mathtext.rm': fam,
        'mathtext.it': f'{fam}:italic',
        'mathtext.bf': f'{fam}:bold',
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.titleweight': 'bold',    # column titles; rows use the y label
        'axes.grid': False,            # pinned: this paper uses a plain ground
        'font.size': 11,
    })


# Colour is keyed to the value of dT over the full sweep, not to its rank in
# whatever subset a figure plots, so a given dT is one colour everywhere.
DT_VMIN, DT_VMAX = 0.2, 1.4
DT_CMAP_FLOOR = 0.0

DT_CMAP = mpl.colors.LinearSegmentedColormap.from_list(
    'viridis_r_trim',
    mpl.colormaps['viridis_r'](np.linspace(DT_CMAP_FLOOR, 1.0, 256)))


def dt_position(dT: float) -> float:
    """Where a divergence sits on the shared colour scale, clipped to the range
        of the standard sweep.
    """
    span = DT_VMAX - DT_VMIN
    return float(np.clip((float(dT) - DT_VMIN) / span, 0.0, 1.0)) if span else 0.5


def dt_colors(task_divs: Sequence[float]) -> Dict[float, tuple]:
    """One colour per task divergence, dark at high divergence."""
    return {dT: DT_CMAP(dt_position(dT)) for dT in task_divs}


def add_dt_colorbar(fig, task_divs: Sequence[float], orientation='horizontal',
                    rect=None):
    norm = mpl.colors.Normalize(vmin=DT_VMIN, vmax=DT_VMAX)
    cmap = DT_CMAP
    if rect is None:
        rect = ([0.27, 0.02, 0.50, 0.020] if orientation == 'horizontal'
                else [0.97, 0.15, 0.015, 0.70])
    ax = fig.add_axes(rect)
    cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm,
                                   orientation=orientation, ticks=task_divs)
    cb.set_ticklabels([f'{dT:.1f}' for dT in task_divs])
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
        ensemble metadata.
    """
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


def wanted_m(m_values, T: int) -> List[int]:
    """Which simultaneities to load for one task number.

        None      every level, 1..T
        'min'     the sequential limit only
        'T'       the simultaneous limit only
        callable  f(T) -> levels
        sequence  explicit levels, dropping any above T
    """
    if m_values is None:
        return list(range(1, T + 1))
    if callable(m_values):
        return sorted({int(m) for m in m_values(T) if m and 1 <= int(m) <= T})
    if isinstance(m_values, str):
        if m_values == 'T':
            return [T]
        if m_values == 'min':
            return [1]
        raise ValueError(f'Unknown m_values selector: {m_values!r}')
    return [int(m) for m in m_values if int(m) <= T]


# Arrays that no metric or figure in this paper reads.
SLIM_DROP = ('P', 'W', 'wait_time', 's_max', 'n_failed_epochs',
             'modularity_entropy', 'cum_time', 'snapshots')

# The complement, passed to `load_condition` so the rest is never decompressed.
SLIM_KEEP = ('pheno_dist', 'd', 'ep_counts', 'active_tasks', 'n_ben')


def slim_reps(reps: List[Dict], drop: Sequence[str] = SLIM_DROP) -> List[Dict]:
    for rep in reps:
        for key in drop:
            rep.pop(key, None)
    return reps


def load_grid(spec: CacheSpec, m_values=None, verbose: bool = True,
              slim: bool = True) -> Dict[int, Dict[float, Dict[int, List]]]:
    """{T: {dT: {m: replicates}}} for everything present in the cache.

        Missing conditions are skipped with a warning; an entirely empty cache
        exits. See `wanted_m` for the selector forms accepted by `m_values`.
    """
    alpha_maps = load_alpha_maps(spec)
    data: Dict[int, Dict[float, Dict[int, List]]] = {}

    for T in spec.T_values:
        if T not in alpha_maps:
            continue
        data[T] = {}
        want_m = wanted_m(m_values, T)

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
                    reps, _ = load_condition(
                        path, keep=SLIM_KEEP if slim else None)
                except (OSError, ValueError, KeyError):
                    continue
                data[T][dT][m] = slim_reps(reps) if slim else reps
            if verbose and data[T][dT]:
                print(f'  T={T} dT={dT}: m={sorted(data[T][dT])} '
                      f'({len(next(iter(data[T][dT].values())))} reps)')

    if not any(ms for ts in data.values() for ms in ts.values()):
        raise SystemExit(
            f'No cached conditions found under {spec.cache_dir}.\n'
            'Build the simulation cache first:  python3 make_data.py')
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
        """Did this replicate reach the requested cutoff, or was it clamped?"""
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
    """Mean pairwise phenotype distance, normalized by this replicate's realized
        task divergence. Bounded in [0, 1].
    """
    td = float(rep['task_dT_realized'])
    if not (np.isfinite(td) and td > 0):
        return np.nan
    return float(rep['pheno_dist'][cutoff.index(rep)]) / td


def optimization(rep: Dict, cutoff: Cutoff) -> float:
    """1 - ||d||_2 / sqrt(T), the fraction of the task deficit eliminated.
        Bounded in [0, 1].
    """
    d = np.asarray(rep['d'], dtype=float)[cutoff.index(rep)]
    return 1.0 - float(np.linalg.norm(d)) / np.sqrt(d.shape[0])


METRICS = {
    'differentiation': (differentiation, r'Degree of differentiation'),
    'optimization': (optimization, r'Degree of optimization'),
}


def metric_values(reps: Sequence[Dict], metric: str,
                  cutoff: Cutoff) -> Tuple[np.ndarray, np.ndarray]:
    """(rep_index, value) arrays for one condition, so that contrasts across
        conditions can be paired on replicate identity.
    """
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
    """Mean and sample SD across replicates."""
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return np.nan, np.nan
    if vals.size == 1:
        return float(vals[0]), 0.0
    return float(vals.mean()), float(vals.std(ddof=1))


def paired_difference(idx_a: np.ndarray, vals_a: np.ndarray,
                      idx_b: np.ndarray, vals_b: np.ndarray
                      ) -> Tuple[float, float, int]:
    """(mean difference, standard error, n) for a - b, matched on rep_index."""
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
    """Termination reasons, realized substitutions, and, if a cutoff is given,
        how many replicates were clamped rather than reaching it.
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
        and the minimum over tasks, both averaged over replicates.
    """
    means, mins = [], []
    for rep in reps:
        ep = np.asarray(rep['ep_counts'], dtype=float)[cutoff.index(rep)]
        means.append(ep.mean())
        mins.append(ep.min())
    if not means:
        return {'ep_mean': np.nan, 'ep_min': np.nan}
    return {'ep_mean': float(np.mean(means)), 'ep_min': float(np.mean(mins))}


def realized_dT_summary(reps: Sequence[Dict]) -> Tuple[float, float]:
    """Mean and SD of the realized task divergence across replicates."""
    return mean_sd(np.array([rep['task_dT_realized'] for rep in reps]))


# ============================================================
# 7. MISC
# ============================================================

def panel_label(i: int) -> str:
    """A, B, ... Z, AA, AB, ..."""
    s = ''
    i += 1
    while i > 0:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


# ============================================================
# 8. VISUAL GRAMMAR
# ============================================================
"""One channel, one meaning, in every main-text and supplementary figure.

  colour      mean pairwise task divergence, dT      (viridis_r, dark = high)
  linestyle   distinguishes series sharing a panel    (dotted 1, dashed T/2,
              -- in practice simultaneity, m            solid T; SOLID whenever
                                                        a panel holds one m)
  marker      number of programs, K                  (only where K varies)
  shaded band dispersion                             (+/- 1 SD, or +/- 1 SE for
                                                      paired differences)

A channel is never used for two things at once, and a quantity already carried
by an axis is not also carried by colour. Reference lines are solid grey, since
dotted now means m = 1.
"""

# --- linestyle: simultaneity -------------------------------------------
LS_M1 = (0, (1, 1.1))        # dense dots; sparse dots vanish at print size
LS_MHALF = (0, (5, 2))
LS_MT = '-'


def ls_for_m(m: Optional[int], T: Optional[int] = None):
    """Linestyle for a simultaneity level. m=None -> solid."""
    if m is None:
        return LS_MT
    if m == 1:
        return LS_M1
    if T is not None and m >= T:
        return LS_MT
    return LS_MHALF


def m_label(m: Optional[int], T: Optional[int] = None) -> str:
    if m is None:
        return ''
    if m == 1:
        return r'$m = 1$'
    if T is not None and m >= T:
        return r'$m = T$'
    return r'$m = T/2$'


# --- marker: number of programs ----------------------------------------
MARKER_FOR_K = {4: 'o', 6: 's', 8: '^'}

# --- dispersion: shaded band -------------------------------------------
def band(ax, x, mu, err, color, ls='-', label=None, lw=1.7, alpha=0.16,
         marker=None, ms=4.5, zorder=2):
    """Mean line over a +/- err band. Bands are drawn below all lines."""
    x = np.asarray(x, dtype=float)
    mu = np.asarray(mu, dtype=float)
    err = np.asarray(err, dtype=float)
    ok = np.isfinite(mu)
    if ok.sum() == 0:
        return None
    ax.fill_between(x[ok], (mu - err)[ok], (mu + err)[ok], color=color,
                    alpha=alpha, lw=0, zorder=zorder)
    line, = ax.plot(x[ok], mu[ok], color=color, ls=ls, lw=lw, label=label,
                    marker=marker, ms=ms if marker else 0,
                    markerfacecolor='none', markeredgecolor=color,
                    markeredgewidth=1.1, zorder=zorder + 10)
    return line


# --- axes ---------------------------------------------------------------
METRIC_YLIM = {'differentiation': (0.0, 1.0), 'optimization': (0.0, 1.0)}


def metric_axis(ax, metric: str, ylabel: bool = True):
    """Fixed limits for a given metric, shared across figures."""
    if metric in METRIC_YLIM:
        ax.set_ylim(*METRIC_YLIM[metric])
    if ylabel:
        ax.set_ylabel(metric_label(metric))


CUTOFF_COLOR = 'darkred'      # comparison points, distinct from the data


def mark_K_line(ax, at: float, label: Optional[str] = None,
                show_label: bool = True, fontsize: float = 11):
    """Solid grey reference at T = K."""
    ax.axvline(at, color='0.68', ls='-', lw=0.9, zorder=0)
    if show_label and label:
        ax.annotate(label, xy=(at, 1.0), xycoords=('data', 'axes fraction'),
                    xytext=(3, -3), textcoords='offset points',
                    fontsize=fontsize, color='0.45', ha='left', va='top')


def mark_cutoff(ax, at: float, label: Optional[str] = None, axis: str = 'x',
                fontsize: float = 10):
    """A comparison point, in dark red so it cannot be read as data."""
    # Above the dispersion bands (zorder 2) but below the mean lines (12), or
    # the band alpha desaturates it to grey.
    (ax.axvline if axis == 'x' else ax.axhline)(
        at, color=CUTOFF_COLOR, ls='-', lw=1.2, zorder=5)
    if not label:
        return
    if axis == 'x':
        ax.annotate(label, xy=(at, 1.0), xycoords=('data', 'axes fraction'),
                    xytext=(3, -3), textcoords='offset points',
                    fontsize=fontsize, color=CUTOFF_COLOR, ha='left', va='top')
    else:
        # right-aligned: the upper left of a panel is usually where the key is
        ax.annotate(label, xy=(0.99, at), xycoords=('axes fraction', 'data'),
                    xytext=(0, 4), textcoords='offset points',
                    fontsize=fontsize, color=CUTOFF_COLOR, ha='right',
                    va='bottom')


def m_legend(ax, levels=('1', 'T/2', 'T'), loc='lower left', fontsize=11,
             **kw):
    """Neutral-grey key for the linestyle channel."""
    from matplotlib.lines import Line2D
    style = {'1': (LS_M1, r'$m = 1$'), 'T/2': (LS_MHALF, r'$m = T/2$'),
             'T': (LS_MT, r'$m = T$')}
    handles = [Line2D([], [], color='0.35', ls=style[k][0], lw=1.7,
                      label=style[k][1]) for k in levels if k in style]
    return ax.legend(handles=handles, loc=loc, fontsize=fontsize,
                     frameon=False, handlelength=2.6, **kw)


