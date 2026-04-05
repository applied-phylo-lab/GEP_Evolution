#!/usr/bin/env python3
"""
fig2_m_sweep.py
===============
Figure 2: endpoint pheno_dist_norm vs m for all m from 1 to T.

Layout: 1 row x 3 cols, one panel per T/K ratio.
x-axis: m (simultaneity, 1 to T)
y-axis: pheno_dist / mean pairwise task-target distance (endpoint mean ± SEM)
Lines:  one per dT value, colored as before (blue=0.2, green=0.6, red=1.2)
Reference line at y=1 (maximum possible differentiation).

Prints numerical table of plotted values.
"""

import json
import os
from itertools import combinations
from typing import Dict, List, Optional

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from simulate import load_condition, load_task_map, sim_cache_path, task_cache_path

mpl.rcParams.update({
    "pdf.use14corefonts": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"],
    "axes.spines.top": False,
    "axes.spines.right": False,
    "font.size": 11,
})

# ── Parameters ───────────────────────────────────────────────

CACHE_DIR = 'simulation_cache'
L, K, GAMMA = 100, 6, 1.0
DENSITY     = 1.0 / K
T_VALUES    = [3, 6, 9]
TASK_DIVS   = [0.2, 0.6, 1.2]
SAVE_DIR    = 'figures_out'
FMT         = 'pdf'

os.makedirs(SAVE_DIR, exist_ok=True)

DT_COLORS = {0.2: '#4878CF', 0.6: '#6ACC65', 1.2: '#D65F5F'}
DT_LABELS = {0.2: r'$\Delta T = 0.2$',
             0.6: r'$\Delta T = 0.6$',
             1.2: r'$\Delta T = 1.2$'}


# ── Metric helpers ────────────────────────────────────────────

def mean_pairwise_distance(X: np.ndarray) -> float:
    """Mean pairwise Euclidean distance between columns of X (L, T)."""
    T = X.shape[1]
    if T < 2:
        return np.nan
    return float(np.mean([
        np.linalg.norm(X[:, i] - X[:, j])
        for i, j in combinations(range(T), 2)
    ]))


# ── Data loading ──────────────────────────────────────────────

def load_all(cache_dir, L, K, gamma, density, T_values, task_divs):
    data, task_maps = {}, {}

    for T in T_values:
        data[T] = {}
        tpath      = task_cache_path(cache_dir, L, K, gamma, T)
        tpath_meta = tpath + '_meta.json'

        task_maps[T] = load_task_map(tpath) \
            if os.path.exists(tpath + '.npz') else {}
        if not task_maps[T]:
            print(f'  WARNING: task map not found for T={T}')

        if not os.path.exists(tpath_meta):
            print(f'  WARNING: task meta not found for T={T}')
            continue

        with open(tpath_meta) as f:
            alpha_map = {float(k): v
                         for k, v in json.load(f)['alpha_map'].items()}

        for dT in task_divs:
            data[T][dT] = {}
            if dT not in alpha_map:
                continue
            alpha = alpha_map[dT]
            for m in range(1, T + 1):
                sp = sim_cache_path(cache_dir, L, K, gamma,
                                    density, T, dT, m, alpha)
                if os.path.exists(sp + '.npz'):
                    results, _ = load_condition(sp)
                    data[T][dT][m] = results

        for dT in task_divs:
            if data[T].get(dT):
                print(f'  T={T} dT={dT}: m={sorted(data[T][dT].keys())}')

    return data, task_maps


# ── Figure ────────────────────────────────────────────────────

def fig2_m_sweep(data: Dict, task_maps: Dict,
                 T_values: List[int], task_divs: List[float], K: int,
                 save_path: Optional[str] = None):
    """
    1 x 3 panel figure. Endpoint pheno_dist_norm vs m.
    Mean ± SEM across replicates. One line per dT.
    """
    n_cols = len(T_values)
    fig, axes = plt.subplots(1, n_cols, figsize=(4.5 * n_cols, 4.5))
    fig.subplots_adjust(wspace=0.15,
                        left=0.10, right=0.96, top=0.88, bottom=0.14)

    # precompute normalization denominators
    target_dists = {
        T: {dT: mean_pairwise_distance(task_maps[T][dT])
            if (T in task_maps and dT in task_maps[T]) else np.nan
            for dT in task_divs}
        for T in T_values
    }

    for col, T in enumerate(T_values):
        ax = axes[col]
        ax.set_box_aspect(1)
        ax.set_title(f'T/K = {T/K:.1f}', fontsize=12, pad=6)

        m_vals_all = sorted(data[T][task_divs[0]].keys()) \
            if task_divs[0] in data[T] else []

        for dT in task_divs:
            if dT not in data[T]:
                continue

            td = target_dists[T][dT]
            if not (np.isfinite(td) and td > 0):
                continue

            means, sems, ms_plot = [], [], []

            for m in sorted(data[T][dT].keys()):
                reps = data[T][dT][m]
                vals = []
                for r in reps:
                    n = r['n_actual_subs']
                    vals.append(float(r['pheno_dist'][n - 1]) / td)
                vals = np.array(vals)
                means.append(vals.mean())
                sems.append(vals.std(ddof=1) / np.sqrt(len(vals)))
                ms_plot.append(m)

            ms_plot = np.array(ms_plot)
            means   = np.array(means)
            sems    = np.array(sems)

            ax.plot(ms_plot, means, '-o',
                    color=DT_COLORS[dT], label=DT_LABELS[dT],
                    lw=1.8, ms=5, zorder=3)
            ax.fill_between(ms_plot,
                            means - sems, means + sems,
                            color=DT_COLORS[dT], alpha=0.15, zorder=2)

        # reference line at maximum possible differentiation
        ax.axhline(1, color='gray', ls='--', lw=0.8, alpha=0.5)

        ax.set_xlabel('m (tasks jointly evaluated)')
        ax.set_ylim(0, None)

        if m_vals_all:
            ax.set_xticks(m_vals_all)

        if col == 0:
            ax.set_ylabel('Realized differentiation\n'
                          r'(pheno dist / task separation)')
        else:
            ax.tick_params(labelleft=False)

    handles, labs = axes[0].get_legend_handles_labels()
    fig.legend(handles, labs,
               loc='upper center', ncol=len(task_divs),
               frameon=False, bbox_to_anchor=(0.5, 0.98), fontsize=10)

    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
        print(f'Saved: {save_path}')

    return fig


# ── Numerical summary ─────────────────────────────────────────

def print_summary(data: Dict, task_maps: Dict,
                  T_values: List[int], task_divs: List[float], K: int):
    """
    Print the exact values plotted: mean ± SEM of endpoint pheno_dist_norm
    for every (T/K, dT, m) combination.
    """
    target_dists = {
        T: {dT: mean_pairwise_distance(task_maps[T][dT])
            if (T in task_maps and dT in task_maps[T]) else np.nan
            for dT in task_divs}
        for T in T_values
    }

    print('\n' + '=' * 60)
    print('NUMERICAL SUMMARY — endpoint pheno_dist_norm')
    print('=' * 60)
    print(f"{'T/K':>5}  {'dT':>5}  {'m':>3}  "
          f"{'mean':>8}  {'sem':>8}  {'n_reps':>6}")
    print('-' * 60)

    for T in T_values:
        for dT in task_divs:
            if dT not in data[T]:
                continue
            td = target_dists[T][dT]
            if not (np.isfinite(td) and td > 0):
                continue
            for m in sorted(data[T][dT].keys()):
                reps = data[T][dT][m]
                vals = np.array([
                    float(r['pheno_dist'][r['n_actual_subs'] - 1]) / td
                    for r in reps
                ])
                print(f"{T/K:>5.1f}  {dT:>5.1f}  {m:>3}  "
                      f"{vals.mean():>8.4f}  "
                      f"{vals.std(ddof=1)/np.sqrt(len(vals)):>8.4f}  "
                      f"{len(vals):>6}")


# ── Run ───────────────────────────────────────────────────────

if __name__ == '__main__':
    print('Loading data...')
    data, task_maps = load_all(
        CACHE_DIR, L, K, GAMMA, DENSITY, T_VALUES, TASK_DIVS)

    out = os.path.join(SAVE_DIR, f'fig2_m_sweep.{FMT}')
    print('Plotting...')
    fig2_m_sweep(data, task_maps, T_VALUES, TASK_DIVS, K, save_path=out)

    print_summary(data, task_maps, T_VALUES, TASK_DIVS, K)

    plt.show()