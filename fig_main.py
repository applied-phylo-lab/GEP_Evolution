#!/usr/bin/env python3
"""
fig_main.py
===========
Main 2x3 figure.

Row 1 — Normalized pheno_dist over log generation time.
         y: pheno_dist / mean pairwise task-target distance
         x: log generation time
         Color family: burgundy = m=1, prussian blue = m=T
         Brightness within family: ΔT=0.2 (light) → ΔT=1.2 (dark)
         Mean trajectory only.

Row 2 — Regulatory architecture start vs end.
         x: mean pairwise Jaccard distance between programs
         y: cosine dissimilarity of program usage
         Start: black asterisk at arrow base (genome identical, usage plastic)
         End: colored dot per (m, ΔT), connected from start by arrow
         Same color scheme as Row 1.

Legend: two colorbars — one burgundy (m=1), one prussian blue (m=T),
        both spanning ΔT = 0.2 → 1.2.

Columns: T/K = 0.5 (T=3), T/K = 1.0 (T=6), T/K = 1.5 (T=9)
"""

import json
import os
from itertools import combinations
from typing import Dict, List, Optional, Tuple

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

# ── Color scheme ─────────────────────────────────────────────
# m=1:  burgundy family
# m=T:  prussian blue family
# ΔT brightness: 0.2=light, 0.6=mid, 1.2=dark (via blend with white)

_BURGUNDY      = np.array([0.549, 0.106, 0.216])
_PRUSSIAN_BLUE = np.array([0.004, 0.267, 0.408])
DT_ALPHAS      = {0.2: 0.30, 0.6: 0.62, 1.2: 1.00}


def condition_color(m_is_T: bool, dT: float) -> tuple:
    base  = _PRUSSIAN_BLUE if m_is_T else _BURGUNDY
    alpha = DT_ALPHAS[dT]
    rgb   = alpha * base + (1 - alpha) * np.ones(3)
    return tuple(rgb)


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


def cosine_dissim(usage: np.ndarray) -> float:
    """Mean pairwise cosine dissimilarity of rows of usage (T, K)."""
    A     = np.array(usage, dtype=float)
    norms = np.linalg.norm(A, axis=1, keepdims=True)
    norms = np.where(norms < 1e-12, 1.0, norms)
    A_n   = A / norms
    sim   = A_n @ A_n.T
    T     = A.shape[0]
    pairs = [(i, j) for i in range(T) for j in range(i + 1, T)]
    if not pairs:
        return 0.0
    return float(1.0 - np.mean([sim[i, j] for i, j in pairs]))


def jaccard_dist(genome: np.ndarray) -> float:
    """Mean pairwise Jaccard distance between program gene sets. genome: (L, K)."""
    K = genome.shape[1]
    if K < 2:
        return 0.0
    dists = []
    for j1, j2 in combinations(range(K), 2):
        g1 = genome[:, j1].astype(bool)
        g2 = genome[:, j2].astype(bool)
        u  = np.sum(g1 | g2)
        dists.append(0.0 if u == 0 else 1.0 - np.sum(g1 & g2) / u)
    return float(np.mean(dists))


def interpolated_mean_logtime(
        reps: List[Dict], y_key: str,
        n_pts: int = 300) -> Tuple[np.ndarray, np.ndarray]:
    """
    Mean trajectory on a common log-generation-time grid.
    Returns (x_grid, mean_y).
    """
    trajs = []
    for r in reps:
        n = r['n_actual_subs']
        x = np.array(r['cum_time'][:n], dtype=float)
        y = np.array(r[y_key][:n],      dtype=float)
        mask = x > 0
        x, y = x[mask], y[mask]
        if len(x) > 1:
            trajs.append((x, y))
    if not trajs:
        return np.array([]), np.array([])

    t_min = min(t[0]  for t, _ in trajs)
    t_max = max(t[-1] for t, _ in trajs)
    if t_min <= 0:
        pos_starts = [t[t > 0][0] for t, _ in trajs if np.any(t > 0)]
        t_min = min(pos_starts) if pos_starts else 1.0
    if t_min >= t_max:
        return np.array([]), np.array([])

    x_grid = np.logspace(np.log10(t_min), np.log10(t_max), n_pts)
    interp  = np.array([np.interp(x_grid, t, v) for t, v in trajs])
    return x_grid, interp.mean(axis=0)


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

def fig_main(data: Dict, task_maps: Dict,
             T_values: List[int], task_divs: List[float], K: int,
             save_path: Optional[str] = None):
    """2 rows x 3 cols main figure."""

    n_cols = len(T_values)
    fig, axes = plt.subplots(2, n_cols, figsize=(4.5 * n_cols, 8.5))
    fig.subplots_adjust(hspace=0.35, wspace=0.20,
                        left=0.10, right=0.92, top=0.88, bottom=0.08)

    # precompute normalization denominators
    target_dists = {
        T: {dT: mean_pairwise_distance(task_maps[T][dT])
                if (T in task_maps and dT in task_maps[T]) else np.nan
            for dT in task_divs}
        for T in T_values
    }

    for col, T in enumerate(T_values):
        ax1 = axes[0, col]
        ax2 = axes[1, col]
        ax1.set_box_aspect(1)
        ax2.set_box_aspect(1)

        ax1.set_title(f'T/K = {T/K:.1f}', fontsize=12, pad=6)

        # ── Compute per-dT initial points ─────────────────────
        # Jaccard of initial genome varies by (rep, m) as noise (~0.008 std)
        # but is statistically exchangeable — average over all m and reps
        # for each dT to get one initial Jaccard per dT.
        # Cosine dissimilarity of initial usage depends on the task geometry
        # (dT), so it genuinely differs across dT conditions. We therefore
        # compute a separate initial point per dT, averaged over all m and reps.
        # Each arrow for a given dT starts from that dT's own initial point.
        # The asterisk(s) mark these per-dT initial points.
        init_by_dt = {}
        for dT in task_divs:
            if dT not in data[T]:
                continue
            j_vals, c_vals = [], []
            for m, reps in data[T][dT].items():
                for r in reps:
                    s0 = sorted(r['snapshots'].keys())[0]
                    g0 = r['snapshots'][s0]['genome']
                    u0 = r['snapshots'][s0]['usage']
                    j_vals.append(jaccard_dist(g0))
                    c_vals.append(cosine_dissim(u0))
            init_by_dt[dT] = (float(np.mean(j_vals)),
                               float(np.mean(c_vals)))

        for dT in task_divs:
            if dT not in data[T]:
                continue

            m_vals = sorted(data[T][dT].keys())
            m1, mT = min(m_vals), max(m_vals)
            td     = target_dists[T][dT]
            ij, ic = init_by_dt[dT]

            for m_is_T, m in [(False, m1), (True, mT)]:
                if m not in data[T][dT]:
                    continue
                reps  = data[T][dT][m]
                color = condition_color(m_is_T, dT)

                # ── Row 1: pheno_dist_norm trajectory ──────────
                # x-axis spans median start to median end time across replicates;
                # replicates with different lengths are silently truncated to
                # this window. Defensible but not neutral — early/late behavior
                # outside the median window is not shown.
                if np.isfinite(td) and td > 0:
                    x_grid, mean_y = interpolated_mean_logtime(
                        reps, y_key='pheno_dist')
                    if len(x_grid) > 0:
                        ax1.plot(x_grid, mean_y / td,
                                 color=color, lw=1.8)

                # ── Row 2: start → end in (Jaccard, cosine) ────
                # Arrow base is the per-dT initial point (ij, ic).
                # Endpoint is mean over replicates of final snapshot.
                e_j, e_c = [], []
                for r in reps:
                    sf = sorted(r['snapshots'].keys())[-1]
                    gf = r['snapshots'][sf]['genome']
                    uf = r['snapshots'][sf]['usage']
                    e_j.append(jaccard_dist(gf))
                    e_c.append(cosine_dissim(uf))

                ex = float(np.mean(e_j))
                ey = float(np.mean(e_c))

                ax2.annotate(
                    '', xy=(ex, ey), xytext=(ij, ic),
                    arrowprops=dict(
                        arrowstyle='->', color=color,
                        lw=1.4, alpha=0.8,
                        mutation_scale=12,
                    )
                )
                ax2.scatter(ex, ey, color=color, s=55, zorder=5,
                            edgecolors='white', linewidths=0.5)

            # asterisk at this dT's initial point
            # color matches the dT brightness using a neutral gray blend
            ax2.scatter(ij, ic, color='black', s=120,
                        zorder=10, marker='*')

        # ── Row 1 axes formatting ──
        ax1.axhline(1, color='gray', ls='--', lw=0.8, alpha=0.5)
        ax1.set_xscale('log')
        ax1.set_ylim(0, None)
        ax1.set_xlabel('Generation (log scale)')
        if col == 0:
            ax1.set_ylabel('Realized differentiation\n'
                           '(pheno dist / task separation)')
        else:
            ax1.tick_params(labelleft=False)

        # ── Row 2 axes formatting ──
        ax2.set_xlim(-0.02, 1.05)
        ax2.set_ylim(-0.02, 1.05)
        ax2.axhline(1, color='gray', ls=':', lw=0.7, alpha=0.35)
        ax2.axvline(1, color='gray', ls=':', lw=0.7, alpha=0.35)
        ax2.set_xlabel('Program individuation\n(mean Jaccard distance)')
        if col == 0:
            ax2.set_ylabel('Usage specialization\n(cosine dissimilarity)')
        else:
            ax2.tick_params(labelleft=False)

    # ── Two colorbars: one per m family ──────────────────────
    # range: 0 to sqrt(2) — natural upper bound for task divergence
    dt_norm = mpl.colors.Normalize(vmin=0, vmax=np.sqrt(2))

    # burgundy colormap (m=1): white → dark burgundy
    cmap_burg = mpl.colors.LinearSegmentedColormap.from_list(
        'burgundy', [np.ones(3), _BURGUNDY])
    # prussian blue colormap (m=T): white → dark prussian blue
    cmap_prus = mpl.colors.LinearSegmentedColormap.from_list(
        'prussian', [np.ones(3), _PRUSSIAN_BLUE])

    # tick positions at the three sampled dT values
    cb_ticks     = task_divs
    cb_ticklabels = [str(dT) for dT in task_divs]

    cbar_y = -0.04
    cbar_h = 0.025
    cbar_w = 0.28

    ax_cb1 = fig.add_axes([0.10, cbar_y, cbar_w, cbar_h])
    ax_cb2 = fig.add_axes([0.55, cbar_y, cbar_w, cbar_h])

    cb1 = mpl.colorbar.ColorbarBase(
        ax_cb1, cmap=cmap_burg, norm=dt_norm,
        orientation='horizontal', ticks=cb_ticks)
    cb1.set_ticklabels(cb_ticklabels)
    cb1.set_label(r'$\Delta T$   (m = 1)', fontsize=9)
    cb1.ax.tick_params(labelsize=8)

    cb2 = mpl.colorbar.ColorbarBase(
        ax_cb2, cmap=cmap_prus, norm=dt_norm,
        orientation='horizontal', ticks=cb_ticks)
    cb2.set_ticklabels(cb_ticklabels)
    cb2.set_label(r'$\Delta T$   (m = T)', fontsize=9)
    cb2.ax.tick_params(labelsize=8)

    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
        print(f'Saved: {save_path}')

    return fig


# ── Run ───────────────────────────────────────────────────────

if __name__ == '__main__':
    print('Loading data...')
    data, task_maps = load_all(
        CACHE_DIR, L, K, GAMMA, DENSITY, T_VALUES, TASK_DIVS)

    out = os.path.join(SAVE_DIR, f'fig_main.{FMT}')
    print('Plotting...')
    fig_main(data, task_maps, T_VALUES, TASK_DIVS, K, save_path=out)

    # ── Numerical summary of what is plotted ──────────────────
    # Row 1: endpoint pheno_dist_norm per (T, dT, m=1, m=T)
    # Row 2: endpoint (Jaccard, cosine) per (T, dT, m=1, m=T)
    #        plus initial (Jaccard, cosine) per (T, dT)
    print('\n' + '='*72)
    print('NUMERICAL SUMMARY')
    print('='*72)

    target_dists = {
        T: {dT: mean_pairwise_distance(task_maps[T][dT])
            if (T in task_maps and dT in task_maps[T]) else np.nan
            for dT in TASK_DIVS}
        for T in T_VALUES
    }

    header = (f"{'T/K':>6}  {'dT':>5}  {'m':>4}  "
              f"{'pheno_norm':>10}  "
              f"{'jaccard_end':>11}  {'cosine_end':>10}  "
              f"{'jaccard_init':>12}  {'cosine_init':>11}")
    print(header)
    print('-' * len(header))

    for T in T_VALUES:
        for dT in TASK_DIVS:
            if dT not in data[T]:
                continue

            td = target_dists[T][dT]

            # initial point for this (T, dT)
            j_init, c_init = [], []
            for m, reps in data[T][dT].items():
                for r in reps:
                    s0 = sorted(r['snapshots'].keys())[0]
                    j_init.append(jaccard_dist(r['snapshots'][s0]['genome']))
                    c_init.append(cosine_dissim(r['snapshots'][s0]['usage']))
            ij = float(np.mean(j_init))
            ic = float(np.mean(c_init))

            m_vals = sorted(data[T][dT].keys())
            m1, mT = min(m_vals), max(m_vals)

            for m in [m1, mT]:
                if m not in data[T][dT]:
                    continue
                reps = data[T][dT][m]

                # pheno_dist_norm at endpoint
                pd_vals = []
                ej_vals, ec_vals = [], []
                for r in reps:
                    n  = r['n_actual_subs']
                    pd_vals.append(float(r['pheno_dist'][n - 1]) / td
                                   if np.isfinite(td) and td > 0 else np.nan)
                    sf = sorted(r['snapshots'].keys())[-1]
                    ej_vals.append(jaccard_dist(r['snapshots'][sf]['genome']))
                    ec_vals.append(cosine_dissim(r['snapshots'][sf]['usage']))

                m_label = '1' if m == m1 else 'T'
                print(f"{T/K:>6.1f}  {dT:>5.1f}  {m_label:>4}  "
                      f"{np.mean(pd_vals):>10.4f}  "
                      f"{np.mean(ej_vals):>11.4f}  {np.mean(ec_vals):>10.4f}  "
                      f"{ij:>12.4f}  {ic:>11.4f}")

    plt.show()