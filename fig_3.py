#!/usr/bin/env python3
"""
fig_3.py
========

pheno_dist_norm vs m sweep at a given substitution cutoff.

Layout: 1 row x N columns, one panel per T/K ratio.
  x-axis: m (simultaneity, 1 to T)
  y-axis: mean pairwise dZ / dT at step `cutoff`
  lines: one per dT, viridis_r color ramp
  ±std stick error bars

Pass --cutoff 50 to inspect the early-time regime.
Default cutoff=None uses the final state (original behaviour).
"""

import argparse
import json
import os
from dataclasses import dataclass, field
from itertools import combinations
from typing import Dict, List, Optional, Tuple

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
    gamma: float = 1.0
    fitness_r: float = 0.0
    densities: List[float] = field(default_factory=lambda: [0.25])
    T_values: List[int] = field(default_factory=lambda: [2, 3, 4, 6, 8])
    task_divs: List[float] = field(default_factory=lambda: [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])
    cutoff: Optional[int] = 200   # None = final state


@dataclass
class OutputConfig:
    save_dir: str = "figures_out"
    filename: str = "fig_3"
    fmt: str = "pdf"
    show: bool = True
    print_summary: bool = True


@dataclass
class FigureConfig:
    figsize_per_col: Tuple[float, float] = (4.5, 4.5)
    line_width: float = 0.75
    marker_size: float = 5.0
    reference_line_y: float = 1.0


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


def pheno_dist_norm_vals(reps: List[Dict], target_dist: float,
                         cutoff: Optional[int]) -> np.ndarray:
    """
    Per-replicate normalized pheno_dist at `cutoff`.
    cutoff=None uses the final state.
    Stalled reps use their actual final value.
    """
    if not (np.isfinite(target_dist) and target_dist > 0):
        return np.array([], dtype=float)
    vals = []
    for r in reps:
        n = r["n_actual_subs"]
        if n <= 0:
            continue
        idx = (min(cutoff - 1, n - 1) if cutoff is not None else n - 1)
        vals.append(float(r["pheno_dist"][idx]) / target_dist)
    return np.array(vals, dtype=float)


def mean_and_std(vals: np.ndarray) -> Tuple[float, float]:
    if len(vals) == 0:
        return np.nan, np.nan
    if len(vals) == 1:
        return float(vals.mean()), 0.0
    return float(vals.mean()), float(vals.std(ddof=1))


def make_dt_viridis_map(task_divs: List[float]) -> Dict[float, Tuple]:
    vals = sorted(task_divs)
    cmap = mpl.colormaps["viridis_r"]
    if len(vals) == 1:
        return {vals[0]: cmap(0.5)}
    return {dT: cmap(i / (len(vals) - 1)) for i, dT in enumerate(vals)}


def _add_viridis_colorbar(fig, task_divs: List[float]):
    dt_norm = mpl.colors.Normalize(vmin=min(task_divs), vmax=max(task_divs))
    cmap = mpl.colormaps["viridis_r"]
    ax_cb = fig.add_axes([0.97, 0.15, 0.015, 0.70])
    cb = mpl.colorbar.ColorbarBase(
        ax_cb, cmap=cmap, norm=dt_norm,
        orientation="vertical", ticks=task_divs
    )
    cb.set_ticklabels([str(dT) for dT in task_divs])
    cb.ax.set_title(r"$\Delta T$", fontsize=10, pad=4)
    cb.ax.tick_params(labelsize=8)


# ============================================================
# 3. DATA LOADING
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
                print(f"WARNING: dT={dT} missing from alpha_map for T={T}")
                continue
            alpha = alpha_map[dT]
            missing_ms = []
            for m in range(1, T + 1):
                sp = sim_cache_path(cfg.cache_dir, cfg.L, cfg.K, cfg.gamma, cfg.fitness_r,
                                    density, T, dT, m, alpha)
                if os.path.exists(sp + ".npz"):
                    results, _ = load_condition(sp)
                    data[T][dT][m] = results
                else:
                    missing_ms.append(m)

            found_ms = sorted(data[T][dT].keys())
            if found_ms:
                print(f"density={density:.4f}  T={T} dT={dT}: m={found_ms}")
            else:
                print(f"WARNING: no data for density={density:.4f}, T={T}, dT={dT}")
            if missing_ms and found_ms:
                print(f"  missing m: {missing_ms}")

    return data, task_maps, alpha_maps


# ============================================================
# 4. FIGURE
# ============================================================

def fig3_m_sweep(
    data: Dict,
    task_maps: Dict,
    data_cfg: DataConfig,
    fig_cfg: FigureConfig,
    save_path: Optional[str] = None,
):
    n_cols = len(data_cfg.T_values)
    fig_w = fig_cfg.figsize_per_col[0] * n_cols
    fig_h = fig_cfg.figsize_per_col[1]
    fig, axes = plt.subplots(1, n_cols, figsize=(fig_w, fig_h))
    if n_cols == 1:
        axes = [axes]

    fig.subplots_adjust(wspace=0.15, left=0.10, right=0.96, top=0.88, bottom=0.14)

    dt_color_map = make_dt_viridis_map(data_cfg.task_divs)
    target_dists = {
        T: {dT: mean_pairwise_distance(task_maps[T][dT])
            if (T in task_maps and dT in task_maps[T]) else np.nan
            for dT in data_cfg.task_divs}
        for T in data_cfg.T_values
    }

    cutoff_label = f"sub ≤ {data_cfg.cutoff}" if data_cfg.cutoff is not None else "final"

    for col, T in enumerate(data_cfg.T_values):
        ax = axes[col]
        ax.set_box_aspect(1)
        ax.set_title(f"T/K = {T / data_cfg.K:.1f}  ({cutoff_label})",
                     fontsize=12, pad=6)

        all_ms = set()

        for dT in data_cfg.task_divs:
            if T not in data or dT not in data[T] or not data[T][dT]:
                continue
            td = target_dists[T][dT]
            if not (np.isfinite(td) and td > 0):
                continue

            ms_plot, means, stds = [], [], []
            for m in sorted(data[T][dT].keys()):
                vals = pheno_dist_norm_vals(data[T][dT][m], td, data_cfg.cutoff)
                if len(vals) == 0:
                    continue
                mean_val, std_val = mean_and_std(vals)
                ms_plot.append(m)
                means.append(mean_val)
                stds.append(std_val)
                all_ms.add(m)

            if not ms_plot:
                continue

            ax.errorbar(
                np.array(ms_plot, dtype=float),
                np.array(means), yerr=np.array(stds),
                fmt="-o", color=dt_color_map[dT],
                lw=fig_cfg.line_width, ms=fig_cfg.marker_size,
                markerfacecolor="none", markeredgecolor=dt_color_map[dT],
                capsize=2, capthick=1.0, elinewidth=1.0, zorder=3,
            )

        ax.axhline(fig_cfg.reference_line_y, color="gray", ls="--", lw=0.8, alpha=0.5)
        ax.set_xlabel("m (tasks jointly evaluated)")
        ax.set_ylim(0, None)
        if all_ms:
            ax.set_xticks(sorted(all_ms))
        if col == 0:
            ax.set_ylabel(r"Mean pairwise $\bar{\Delta Z} / \bar{\Delta T}$")
        else:
            ax.tick_params(labelleft=False)

    _add_viridis_colorbar(fig, data_cfg.task_divs)

    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
        print(f"Saved: {save_path}")
    return fig


# ============================================================
# 5. NUMERICAL SUMMARY
# ============================================================

def print_summary(data: Dict, task_maps: Dict, data_cfg: DataConfig):
    target_dists = {
        T: {dT: mean_pairwise_distance(task_maps[T][dT])
            if (T in task_maps and dT in task_maps[T]) else np.nan
            for dT in data_cfg.task_divs}
        for T in data_cfg.T_values
    }

    cutoff_label = f"sub ≤ {data_cfg.cutoff}" if data_cfg.cutoff is not None else "final"
    print(f"\n{'=' * 68}")
    print(f"NUMERICAL SUMMARY — pheno_dist_norm ({cutoff_label})")
    print("=" * 68)
    print(f"{'T/K':>5}  {'dT':>5}  {'m':>3}  {'mean':>8}  {'std':>8}  {'n_reps':>6}")
    print("-" * 68)

    for T in data_cfg.T_values:
        for dT in data_cfg.task_divs:
            if T not in data or dT not in data[T] or not data[T][dT]:
                continue
            td = target_dists[T][dT]
            if not (np.isfinite(td) and td > 0):
                continue
            for m in sorted(data[T][dT].keys()):
                vals = pheno_dist_norm_vals(data[T][dT][m], td, data_cfg.cutoff)
                if len(vals) == 0:
                    continue
                mean_val, std_val = mean_and_std(vals)
                print(f"{T / data_cfg.K:>5.1f}  {dT:>5.1f}  {m:>3}  "
                      f"{mean_val:>8.4f}  {std_val:>8.4f}  {len(vals):>6}")


# ============================================================
# 6. CLI + MAIN
# ============================================================

def parse_args():
    d = DataConfig()
    o = OutputConfig()
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
    p.add_argument("--cutoff",     type=int,   default=d.cutoff)
    p.add_argument("--no_show",    action="store_true")
    p.add_argument("--no_summary", action="store_true")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()

    data_cfg = DataConfig(
        cache_dir=args.cache_dir, L=args.L, K=args.K, gamma=args.gamma,
        fitness_r=args.fitness_r,
        densities=args.densities, T_values=args.T_values,
        task_divs=args.task_divs, cutoff=args.cutoff,
    )
    out_cfg = OutputConfig(
        save_dir=args.save_dir, filename=args.filename,
        fmt=args.fmt, show=not args.no_show,
        print_summary=not args.no_summary,
    )
    fig_cfg = FigureConfig()
    os.makedirs(out_cfg.save_dir, exist_ok=True)

    for density in data_cfg.densities:
        print(f"\nLoading data for density={density:.4f}...")
        data, task_maps, _ = load_all_for_density(data_cfg, density)

        cutoff_tag = f"_sub{data_cfg.cutoff}" if data_cfg.cutoff is not None else ""
        save_path = os.path.join(
            out_cfg.save_dir,
            f"{out_cfg.filename}{cutoff_tag}"
            f"_gamma{data_cfg.gamma}_fr{data_cfg.fitness_r}"
            f"_density{density:.4f}.{out_cfg.fmt}"
        )
        print(f"Plotting density={density:.4f}...")
        fig = fig3_m_sweep(data, task_maps, data_cfg, fig_cfg, save_path=save_path)

        if out_cfg.print_summary:
            print_summary(data, task_maps, data_cfg)

        if out_cfg.show:
            plt.show()
        else:
            plt.close(fig)