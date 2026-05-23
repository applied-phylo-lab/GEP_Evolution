#!/usr/bin/env python3
"""
fig_s_shape_inheritance.py
==========================

Supplementary figure testing the claim:

    "Intermediate m shows T/K dependence more like m=T than like m=1."

The two reference panels (Figure 2 of the main text) plot Delta_Z / Delta_T
versus T/K, with one line per Delta_T, in two panels:
    Panel A (m = 1):     monotone decline across T/K, including T/K < 1.
    Panel B (m = T):     plateau for T/K <= 1, decline only for T/K > 1.

This figure replots the same quantity with one line per *normalized simultaneity*
m_norm = (m - 1) / (T - 1), so lines are directly comparable across T values:
    m_norm = 0    <-> m = 1     (strictly sequential)
    m_norm = 1    <-> m = T     (fully simultaneous)
    intermediate  <-> partial simultaneity

m_norm is binned into a small number of bands (default 4) and lines show the
mean over m values that fall in each band. One panel per Delta_T.

The claim is supported if intermediate-band lines sit visibly closer to the
m_norm = 1 line than to the m_norm = 0 line in the T/K < 1 regime, i.e. if
the plateau-then-decline shape appears as soon as m_norm > 0 rather than only
at m_norm = 1.

USAGE
    python fig_s_shape_inheritance.py
    python fig_s_shape_inheritance.py --cutoffs 100 200 --n_bands 5
    python fig_s_shape_inheritance.py --dT 0.4 0.8 1.2  # subset of Delta_T panels
"""

import argparse
import json
import os
from dataclasses import dataclass, field
from itertools import combinations
from typing import Dict, List, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from simulate import load_condition, load_task_map, sim_cache_path, task_cache_path


# ============================================================
# 0. MPL STYLE
# ============================================================

mpl.rcParams.update({
    "pdf.use14corefonts": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"],
    "axes.spines.top": False,
    "axes.spines.right": False,
    "font.size": 11,
})


# ============================================================
# 1. CONFIG
# ============================================================

@dataclass
class DataConfig:
    cache_dir: str = "simulation_cache"
    L: int = 100
    K: int = 4
    gamma: float = 4.0
    fitness_r: float = -2
    densities: List[float] = field(default_factory=lambda: [1/4])
    T_values: List[int] = field(default_factory=lambda: [2, 3, 4, 6, 8])
    task_divs: List[float] = field(default_factory=lambda: [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])
    cutoffs: List[int] = field(default_factory=lambda: [200])


@dataclass
class OutputConfig:
    save_dir: str = "figures_out"
    filename: str = "fig_s_shape_inheritance"
    fmt: str = "pdf"
    show: bool = True
    print_summary: bool = True


@dataclass
class FigureConfig:
    # one panel per Delta_T value (subset selectable on command line)
    panel_task_divs: List[float] = field(default_factory=lambda: [1.2,1.4])
    n_bands: int = 4
    panel_size: Tuple[float, float] = (3.2, 3.2)


# ============================================================
# 2. HELPERS
# ============================================================

def mean_pairwise_distance(X: np.ndarray) -> float:
    T = X.shape[1]
    if T < 2:
        return np.nan
    return float(np.mean([
        np.linalg.norm(X[:, i] - X[:, j])
        for i, j in combinations(range(T), 2)
    ]))


def pheno_dist_norm_at_step(reps: List[Dict], target_dist: float, cutoff: int) -> np.ndarray:
    if not np.isfinite(target_dist) or target_dist <= 0:
        return np.array([])
    return np.array([
        float(r["pheno_dist"][min(cutoff - 1, r["n_actual_subs"] - 1)]) / target_dist
        for r in reps
    ])


def mean_and_std(vals: np.ndarray) -> Tuple[float, float]:
    if len(vals) == 0:
        return np.nan, np.nan
    if len(vals) == 1:
        return float(vals.mean()), 0.0
    return float(vals.mean()), float(vals.std(ddof=1))


def m_norm(m: int, T: int) -> float:
    """Normalized simultaneity: 0 for m=1, 1 for m=T."""
    if T <= 1:
        return 0.0
    return (m - 1) / (T - 1)


def band_of(m_norm_val: float, n_bands: int) -> int:
    """
    Assign m_norm in [0, 1] to one of n_bands bands.

    Endpoints are pinned: m_norm == 0 -> band 0, m_norm == 1 -> band n_bands-1.
    Intermediate values fall into evenly spaced interior bands.
    """
    if m_norm_val <= 0.0:
        return 0
    if m_norm_val >= 1.0:
        return n_bands - 1
    # interior values into bands 0..n_bands-1 by simple binning
    idx = int(np.floor(m_norm_val * n_bands))
    return min(idx, n_bands - 1)


def band_label(b: int, n_bands: int) -> str:
    """Pretty label for a band of normalized simultaneity."""
    if b == 0:
        return r"$m{=}1$ (seq.)"
    if b == n_bands - 1:
        return r"$m{=}T$ (sim.)"
    lo = b / n_bands
    hi = (b + 1) / n_bands
    return f"$m_{{norm}} \\in [{lo:.2f}, {hi:.2f}]$"


def make_band_cmap(n_bands: int):
    cmap = mpl.colormaps["plasma"]
    return [cmap(i / max(n_bands - 1, 1)) for i in range(n_bands)]


# ============================================================
# 3. DATA LOADING (identical to fig_2.py)
# ============================================================

def load_all_for_density(cfg: DataConfig, density: float):
    data, task_maps, alpha_maps = {}, {}, {}

    for T in cfg.T_values:
        data[T] = {}
        tpath = task_cache_path(cfg.cache_dir, cfg.L, cfg.K, cfg.gamma, cfg.fitness_r, T)
        tpath_meta = tpath + "_meta.json"

        task_maps[T] = load_task_map(tpath) if os.path.exists(tpath + ".npz") else {}
        if not task_maps[T]:
            print(f"WARNING: task map not found for T={T}")
        if not os.path.exists(tpath_meta):
            print(f"WARNING: task meta not found for T={T}")
            continue

        with open(tpath_meta) as f:
            alpha_map = {float(k): v for k, v in json.load(f)["alpha_map"].items()}
        alpha_maps[T] = alpha_map

        for dT in cfg.task_divs:
            data[T][dT] = {}
            if dT not in alpha_map:
                continue
            alpha = alpha_map[dT]
            for m in range(1, T + 1):
                sp = sim_cache_path(cfg.cache_dir, cfg.L, cfg.K, cfg.gamma, cfg.fitness_r,
                                    density, T, dT, m, alpha)
                if os.path.exists(sp + ".npz"):
                    results, _ = load_condition(sp)
                    data[T][dT][m] = results

    return data, task_maps, alpha_maps


# ============================================================
# 4. BAND AGGREGATION
# ============================================================

def aggregate_by_band(
    data, task_maps, data_cfg, fig_cfg, cutoff, dT
):
    """
    For a given Delta_T value, return:
        x_by_band:    {band -> array of T/K values}
        y_mean_by_band: {band -> array of mean Delta_Z/Delta_T}
        y_std_by_band:  {band -> array of std across replicates pooled within band}

    For each (T, m) cell, compute replicate-level normalized differentiation.
    Then for each T, pool replicates from all m values assigned to the same band
    and report the band's (mean, std) at that T.
    """
    n_bands = fig_cfg.n_bands

    x_by_band = {b: [] for b in range(n_bands)}
    y_mean_by_band = {b: [] for b in range(n_bands)}
    y_std_by_band = {b: [] for b in range(n_bands)}

    for T in data_cfg.T_values:
        if dT not in data.get(T, {}) or not data[T][dT]:
            continue
        if T not in task_maps or dT not in task_maps[T]:
            continue
        target_dist = mean_pairwise_distance(task_maps[T][dT])

        # pool replicate values across all m in each band
        band_vals = {b: [] for b in range(n_bands)}
        for m, reps in data[T][dT].items():
            mn = m_norm(m, T)
            b = band_of(mn, n_bands)
            vals = pheno_dist_norm_at_step(reps, target_dist, cutoff)
            if len(vals) > 0:
                band_vals[b].extend(vals.tolist())

        tk = T / data_cfg.K
        for b in range(n_bands):
            if band_vals[b]:
                arr = np.array(band_vals[b])
                m_, s_ = mean_and_std(arr)
                x_by_band[b].append(tk)
                y_mean_by_band[b].append(m_)
                y_std_by_band[b].append(s_)

    # to arrays
    for b in range(n_bands):
        x_by_band[b] = np.array(x_by_band[b])
        y_mean_by_band[b] = np.array(y_mean_by_band[b])
        y_std_by_band[b] = np.array(y_std_by_band[b])

    return x_by_band, y_mean_by_band, y_std_by_band


# ============================================================
# 5. FIGURE
# ============================================================

def fig_shape_inheritance(data, task_maps, data_cfg, fig_cfg, cutoff, save_path=None):
    """One panel per Delta_T: Delta_Z/Delta_T vs T/K, lines indexed by simultaneity band."""
    panel_dTs = [d for d in fig_cfg.panel_task_divs if d in data_cfg.task_divs]
    if not panel_dTs:
        raise ValueError("None of the requested panel Delta_T values are in task_divs.")

    n_panels = len(panel_dTs)
    fig, axes = plt.subplots(
        1, n_panels,
        figsize=(fig_cfg.panel_size[0] * n_panels + 1.0, fig_cfg.panel_size[1] + 1.0),
        sharey=True,
    )
    if n_panels == 1:
        axes = np.array([axes])

    fig.subplots_adjust(wspace=0.18, left=0.10, right=0.78, top=0.86, bottom=0.18)

    n_bands = fig_cfg.n_bands
    band_colors = make_band_cmap(n_bands)
    panel_labels = "ABCDEFGH"

    tk_values = np.array(sorted({T / data_cfg.K for T in data_cfg.T_values}))

    for pi, dT in enumerate(panel_dTs):
        ax = axes[pi]
        ax.set_box_aspect(1)
        ax.text(-0.18, 1.08, panel_labels[pi], transform=ax.transAxes,
                fontsize=14, fontweight="bold", va="top", ha="left")
        ax.set_title(rf"$\Delta T = {dT}$  (sub $\leq$ {cutoff})", fontsize=11, pad=6)

        x_by_band, ym_by_band, ys_by_band = aggregate_by_band(
            data, task_maps, data_cfg, fig_cfg, cutoff, dT
        )

        for b in range(n_bands):
            if len(x_by_band[b]) == 0:
                continue
            # sort by T/K so the line is monotone in x
            order = np.argsort(x_by_band[b])
            x = x_by_band[b][order]
            ym = ym_by_band[b][order]
            ys = ys_by_band[b][order]

            ax.errorbar(
                x, ym, yerr=ys,
                fmt="-o", color=band_colors[b], lw=1.0, ms=5,
                markerfacecolor="none", markeredgecolor=band_colors[b],
                capsize=2, capthick=1.0, elinewidth=1.0,
                label=band_label(b, n_bands),
            )

        ax.axhline(1.0, color="gray", ls="--", lw=0.8, alpha=0.5)
        ax.axvline(1.0, color="gray", ls=":", lw=0.8, alpha=0.4)
        ax.set_xticks(tk_values)
        ax.set_xticklabels([f"{v:.1f}" for v in tk_values])
        ax.set_xlabel("T / K")
        ax.set_ylim(0, None)
        if pi == 0:
            ax.set_ylabel(r"Mean pairwise $\bar{\Delta Z} / \bar{\Delta T}$")

    # legend outside, right of all panels
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(
        handles, labels,
        loc="center left",
        bbox_to_anchor=(0.80, 0.5),
        frameon=False,
        fontsize=9,
        title="simultaneity band",
        title_fontsize=9,
    )

    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
        print(f"Saved: {save_path}")
    return fig


# ============================================================
# 6. NUMERICAL SUMMARY
# ============================================================

def print_shape_summary(data, task_maps, data_cfg, fig_cfg, cutoff):
    """
    For each Delta_T in panels, print the gap between the m=1 band and each
    higher band, at each T/K. Useful for showing quantitatively that
    intermediate-band curves lift off the m=1 curve immediately.
    """
    panel_dTs = [d for d in fig_cfg.panel_task_divs if d in data_cfg.task_divs]
    n_bands = fig_cfg.n_bands

    print(f"\n{'=' * 96}")
    print(f"SHAPE-INHERITANCE SUMMARY — sub ≤ {cutoff}")
    print("(values are mean Delta_Z / Delta_T; gap = band b minus band 0 at same T/K)")
    print("=" * 96)

    for dT in panel_dTs:
        print(f"\nDelta_T = {dT}")
        x_by_band, ym_by_band, _ = aggregate_by_band(
            data, task_maps, data_cfg, fig_cfg, cutoff, dT
        )

        # build a tk-indexed lookup per band
        tk_set = sorted({tk for b in range(n_bands) for tk in x_by_band[b].tolist()})
        lookup = {b: dict(zip(x_by_band[b].tolist(), ym_by_band[b].tolist()))
                  for b in range(n_bands)}

        header = f"  {'T/K':>5}  " + "  ".join(
            f"{band_label(b, n_bands).replace('$', '').replace('{', '').replace('}', ''):>22}"
            for b in range(n_bands)
        )
        print(header)
        for tk in tk_set:
            row = [f"  {tk:>5.2f}"]
            base = lookup[0].get(tk, None)
            for b in range(n_bands):
                v = lookup[b].get(tk, None)
                if v is None:
                    row.append(f"{'--':>22}")
                else:
                    if b == 0 or base is None:
                        row.append(f"{v:>10.3f}  (gap  N/A)")
                    else:
                        gap = v - base
                        row.append(f"{v:>10.3f}  (gap {gap:>+5.3f})")
            print("  ".join(row))


# ============================================================
# 7. CLI + MAIN
# ============================================================

def parse_args():
    d = DataConfig()
    o = OutputConfig()
    f = FigureConfig()
    p = argparse.ArgumentParser()
    p.add_argument("--cache_dir",  type=str,   default=d.cache_dir)
    p.add_argument("--save_dir",   type=str,   default=o.save_dir)
    p.add_argument("--filename",   type=str,   default=o.filename)
    p.add_argument("--fmt",        type=str,   default=o.fmt)
    p.add_argument("--L",          type=int,   default=d.L)
    p.add_argument("--K",          type=int,   default=d.K)
    p.add_argument("--gamma",      type=float, default=d.gamma)
    p.add_argument("--fitness_r",  type=float, default=d.fitness_r)
    p.add_argument("--densities",  type=float, nargs="+", default=d.densities)
    p.add_argument("--T",          type=int,   nargs="+", default=d.T_values,  dest="T_values")
    p.add_argument("--dT",         type=float, nargs="+", default=d.task_divs, dest="task_divs")
    p.add_argument("--panel_dT",   type=float, nargs="+", default=f.panel_task_divs)
    p.add_argument("--n_bands",    type=int,   default=f.n_bands)
    p.add_argument("--cutoffs",    type=int,   nargs="+", default=d.cutoffs)
    p.add_argument("--no_show",    action="store_true")
    p.add_argument("--no_summary", action="store_true")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()

    data_cfg = DataConfig(
        cache_dir=args.cache_dir, L=args.L, K=args.K, gamma=args.gamma,
        fitness_r=args.fitness_r,
        densities=args.densities, T_values=args.T_values,
        task_divs=args.task_divs, cutoffs=args.cutoffs,
    )
    out_cfg = OutputConfig(
        save_dir=args.save_dir, filename=args.filename,
        fmt=args.fmt, show=not args.no_show,
        print_summary=not args.no_summary,
    )
    fig_cfg = FigureConfig(
        panel_task_divs=args.panel_dT,
        n_bands=args.n_bands,
    )
    os.makedirs(out_cfg.save_dir, exist_ok=True)

    for density in data_cfg.densities:
        print(f"\nLoading data for density={density:.4f}...")
        data, task_maps, _ = load_all_for_density(data_cfg, density)

        for cutoff in data_cfg.cutoffs:
            def sp(tag):
                return os.path.join(
                    out_cfg.save_dir,
                    f"{out_cfg.filename}_{tag}_sub{cutoff}"
                    f"_gamma{data_cfg.gamma}_fr{data_cfg.fitness_r}"
                    f"_density{density:.4f}_bands{fig_cfg.n_bands}.{out_cfg.fmt}"
                )

            fig = fig_shape_inheritance(
                data, task_maps, data_cfg, fig_cfg,
                cutoff=cutoff, save_path=sp("main"),
            )

            if out_cfg.print_summary:
                print_shape_summary(data, task_maps, data_cfg, fig_cfg, cutoff)

            if out_cfg.show:
                plt.show()
            else:
                plt.close(fig)