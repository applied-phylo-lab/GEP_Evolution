#!/usr/bin/env python3
"""
fig_trajectories.py
===================

Substitution trajectories of mean pairwise dZ/dT.

Figure 1: N_T rows x 2 columns
  Columns: m=1 (left), m=T (right)
  Rows: one per T/K ratio
  Each panel: all individual replicate trajectories, colored by dT (viridis_r).
  Horizontal colorbar at bottom.

Figure 2a/2b: 1 row x N_m columns, one figure per target T/K
  Columns: all available m values
  Each panel: all replicate trajectories colored by dT.
  Vertical colorbar on the right.

x-axis: substitution count (downsampled 1/10)
y-axis: mean pairwise dZ/dT
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


@dataclass
class OutputConfig:
    save_dir: str = "figures_out"
    filename: str = "fig_trajectories"
    fmt: str = "pdf"
    show: bool = True


@dataclass
class FigureConfig:
    figsize_per_panel: Tuple[float, float] = (3.5, 3.0)
    line_width: float = 0.4
    alpha: float = 0.5


# ============================================================
# 2. METRIC HELPERS
# ============================================================

def mean_pairwise_distance(X: np.ndarray) -> float:
    T = X.shape[1]
    if T < 2:
        return np.nan
    return float(np.mean([
        np.linalg.norm(X[:, i] - X[:, j])
        for i, j in combinations(range(T), 2)
    ]))


# ============================================================
# 3. COLOR / COLORBAR HELPERS
# ============================================================

def make_dt_viridis_map(task_divs: List[float]) -> Dict[float, Tuple]:
    """Map each dT to viridis_r color, linearly spaced."""
    vals = sorted(task_divs)
    cmap = mpl.colormaps["viridis_r"]
    if len(vals) == 1:
        return {vals[0]: cmap(0.5)}
    return {dT: cmap(i / (len(vals) - 1)) for i, dT in enumerate(vals)}


def _add_colorbar_horizontal(fig, task_divs: List[float]):
    """Horizontal viridis_r colorbar centered at the bottom."""
    dt_norm = mpl.colors.Normalize(vmin=min(task_divs), vmax=max(task_divs))
    cmap = mpl.colormaps["viridis_r"]
    cb_w = 0.50
    ax_cb = fig.add_axes([0.5 - cb_w / 2 + 0.02, 0.02, cb_w, 0.018])
    cb = mpl.colorbar.ColorbarBase(
        ax_cb, cmap=cmap, norm=dt_norm,
        orientation="horizontal", ticks=task_divs
    )
    cb.set_ticklabels([str(dT) for dT in task_divs])
    cb.set_label(r"$\Delta T$", fontsize=10)
    cb.ax.tick_params(labelsize=8)


def _add_colorbar_vertical(fig, task_divs: List[float]):
    """Vertical viridis_r colorbar on the right, Delta T label on top."""
    dt_norm = mpl.colors.Normalize(vmin=min(task_divs), vmax=max(task_divs))
    cmap = mpl.colormaps["viridis_r"]
    ax_cb = fig.add_axes([0.95, 0.15, 0.015, 0.70])
    cb = mpl.colorbar.ColorbarBase(
        ax_cb, cmap=cmap, norm=dt_norm,
        orientation="vertical", ticks=task_divs
    )
    cb.set_ticklabels([str(dT) for dT in task_divs])
    cb.ax.set_title(r"$\Delta T$", fontsize=10, pad=4)
    cb.ax.tick_params(labelsize=8)


# ============================================================
# 4. DATA LOADING
# ============================================================

def load_all_for_density(cfg: DataConfig, density: float):
    data = {}
    task_maps = {}
    alpha_maps = {}

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
                sp = sim_cache_path(
                    cfg.cache_dir, cfg.L, cfg.K, cfg.gamma, cfg.fitness_r,
                    density, T, dT, m, alpha
                )
                if os.path.exists(sp + ".npz"):
                    results, _ = load_condition(sp)
                    data[T][dT][m] = results

        for dT in cfg.task_divs:
            if data[T].get(dT):
                print(f"density={density:.4f}  T={T} dT={dT}: m={sorted(data[T][dT].keys())}")

    return data, task_maps, alpha_maps


# ============================================================
# 5. SHARED PANEL PLOTTER
# ============================================================

def _plot_panel(ax, reps_by_dT, td_map, dt_color_map, fig_cfg):
    """Plot all replicate trajectories (downsampled 1/10) into a single axis."""
    for dT, reps in reps_by_dT.items():
        td = td_map.get(dT, np.nan)
        if not (np.isfinite(td) and td > 0):
            continue
        color = dt_color_map[dT]
        for r in reps:
            n = r["n_actual_subs"]
            x = np.arange(n)[::10]
            y = np.array(r["pheno_dist"][:n], dtype=float)[::10] / td
            ax.plot(x, y, color=color, lw=fig_cfg.line_width, alpha=fig_cfg.alpha)
    ax.axhline(1, color="gray", ls="--", lw=0.8, alpha=0.5)
    ax.set_ylim(0, None)


# ============================================================
# 6. FIGURE 1 — rows=T/K, cols=m=1|m=T, horizontal colorbar
# ============================================================

def fig_trajectories(
    data: Dict,
    task_maps: Dict,
    data_cfg: DataConfig,
    fig_cfg: FigureConfig,
    save_path: Optional[str] = None,
):
    """N_T rows x 2 cols. Cols = m=1, m=T. Rows = T/K. Horizontal colorbar."""
    n_rows = len(data_cfg.T_values)

    fig_w = fig_cfg.figsize_per_panel[0] * 2
    fig_h = fig_cfg.figsize_per_panel[1] * n_rows

    fig, axes = plt.subplots(n_rows, 2, figsize=(fig_w, fig_h))
    if n_rows == 1:
        axes = axes[np.newaxis, :]

    fig.subplots_adjust(
        wspace=0.10, hspace=0.30,
        left=0.12, right=0.96, top=0.95, bottom=0.08
    )

    dt_color_map = make_dt_viridis_map(data_cfg.task_divs)
    target_dists = {
        T: {
            dT: mean_pairwise_distance(task_maps[T][dT])
            if (T in task_maps and dT in task_maps[T]) else np.nan
            for dT in data_cfg.task_divs
        }
        for T in data_cfg.T_values
    }

    panel_labels = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    pi = 0

    for row, T in enumerate(data_cfg.T_values):
        td_map = target_dists[T]

        for col, (m, col_label) in enumerate([(1, "m = 1"), (T, "m = T")]):
            ax = axes[row, col]
            ax.set_box_aspect(1)

            ax.text(-0.15, 1.08, panel_labels[pi],
                    transform=ax.transAxes, fontsize=14, fontweight="bold",
                    va="top", ha="left")
            pi += 1

            if row == 0:
                ax.set_title(col_label, fontsize=12, pad=6)

            if col == 0:
                ax.set_ylabel(
                    f"T/K = {T / data_cfg.K:.1f}\n"
                    + r"$\bar{\Delta Z} / \bar{\Delta T}$",
                    fontsize=10,
                )
            else:
                ax.tick_params(labelleft=False)

            if row == n_rows - 1:
                ax.set_xlabel("Substitution count")
            else:
                ax.tick_params(labelbottom=False)

            reps_by_dT = {
                dT: data[T][dT][m]
                for dT in data_cfg.task_divs
                if dT in data[T] and m in data[T][dT]
            }
            _plot_panel(ax, reps_by_dT, td_map, dt_color_map, fig_cfg)

    _add_colorbar_horizontal(fig, data_cfg.task_divs)

    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
        print(f"Saved: {save_path}")

    return fig


# ============================================================
# 7. FIGURE 2 — cols=m, single T/K, vertical colorbar
# ============================================================

def fig_trajectories_m_sweep(
    data: Dict,
    task_maps: Dict,
    data_cfg: DataConfig,
    fig_cfg: FigureConfig,
    target_tk: float = 0.8,
    save_path: Optional[str] = None,
):
    """1 row x N_m cols for the T/K closest to target_tk. Vertical colorbar."""
    T = min(data_cfg.T_values, key=lambda t: abs(t / data_cfg.K - target_tk))
    actual_tk = T / data_cfg.K

    all_m = sorted({
        m for dT in data_cfg.task_divs
        if dT in data[T]
        for m in data[T][dT].keys()
    })
    if not all_m:
        print(f"No data found for T={T}")
        return None

    n_cols = len(all_m)
    fig_w = fig_cfg.figsize_per_panel[0] * n_cols
    fig_h = fig_cfg.figsize_per_panel[1]

    fig, axes = plt.subplots(1, n_cols, figsize=(fig_w, fig_h))
    if n_cols == 1:
        axes = [axes]

    fig.subplots_adjust(
        wspace=0.10,
        left=0.12, right=0.93, top=0.88, bottom=0.12
    )

    dt_color_map = make_dt_viridis_map(data_cfg.task_divs)
    td_map = {
        dT: mean_pairwise_distance(task_maps[T][dT])
        if (T in task_maps and dT in task_maps[T]) else np.nan
        for dT in data_cfg.task_divs
    }
    panel_labels = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    for col, m in enumerate(all_m):
        ax = axes[col]
        ax.set_box_aspect(1)

        ax.text(-0.15, 1.08, panel_labels[col],
                transform=ax.transAxes, fontsize=14, fontweight="bold",
                va="top", ha="left")
        ax.set_title(f"m = {m}", fontsize=12, pad=6)
        ax.set_xlabel("Substitution count")

        if col == 0:
            ax.set_ylabel(r"$\bar{\Delta Z} / \bar{\Delta T}$", fontsize=10)
        else:
            ax.tick_params(labelleft=False)

        reps_by_dT = {
            dT: data[T][dT][m]
            for dT in data_cfg.task_divs
            if dT in data[T] and m in data[T][dT]
        }
        _plot_panel(ax, reps_by_dT, td_map, dt_color_map, fig_cfg)

    fig.suptitle(f"T/K = {actual_tk:.1f}", fontsize=13, y=1.02)
    _add_colorbar_vertical(fig, data_cfg.task_divs)

    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
        print(f"Saved: {save_path}")

    return fig


# ============================================================
# 8. CLI
# ============================================================

def parse_args():
    d_cfg = DataConfig()
    o_cfg = OutputConfig()
    p = argparse.ArgumentParser()

    p.add_argument("--cache_dir",  type=str,   default=d_cfg.cache_dir)
    p.add_argument("--save_dir",   type=str,   default=o_cfg.save_dir)
    p.add_argument("--filename",   type=str,   default=o_cfg.filename)
    p.add_argument("--fmt",        type=str,   default=o_cfg.fmt)
    p.add_argument("--L",          type=int,   default=d_cfg.L)
    p.add_argument("--K",          type=int,   default=d_cfg.K)
    p.add_argument("--gamma",      type=float, default=d_cfg.gamma)
    p.add_argument("--fitness_r",  type=float, default=d_cfg.fitness_r)
    p.add_argument("--densities",  type=float, nargs="+", default=d_cfg.densities)
    p.add_argument("--T",          type=int,   nargs="+", default=d_cfg.T_values,  dest="T_values")
    p.add_argument("--dT",         type=float, nargs="+", default=d_cfg.task_divs, dest="task_divs")
    p.add_argument("--no_show",    action="store_true")

    return p.parse_args()


# ============================================================
# 9. MAIN
# ============================================================

if __name__ == "__main__":
    args = parse_args()

    data_cfg = DataConfig(
        cache_dir=args.cache_dir,
        L=args.L,
        K=args.K,
        gamma=args.gamma,
        fitness_r=args.fitness_r,
        densities=args.densities,
        T_values=args.T_values,
        task_divs=args.task_divs,
    )

    out_cfg = OutputConfig(
        save_dir=args.save_dir,
        filename=args.filename,
        fmt=args.fmt,
        show=not args.no_show,
    )

    fig_cfg = FigureConfig()
    os.makedirs(out_cfg.save_dir, exist_ok=True)

    for density in data_cfg.densities:
        print(f"\nLoading data for density={density:.4f}...")
        data, task_maps, _ = load_all_for_density(data_cfg, density)

        tag = f"_gamma{data_cfg.gamma}_fr{data_cfg.fitness_r}_density{density:.4f}"

        # Figure 1: rows=T/K, cols=m=1|m=T
        fig = fig_trajectories(
            data=data, task_maps=task_maps,
            data_cfg=data_cfg, fig_cfg=fig_cfg,
            save_path=os.path.join(
                out_cfg.save_dir,
                f"{out_cfg.filename}{tag}.{out_cfg.fmt}"
            ),
        )

        # Figure 2a: m sweep at T/K = 0.8
        fig_m1 = fig_trajectories_m_sweep(
            data=data, task_maps=task_maps,
            data_cfg=data_cfg, fig_cfg=fig_cfg,
            target_tk=0.8,
            save_path=os.path.join(
                out_cfg.save_dir,
                f"{out_cfg.filename}_m_sweep_tk0.8{tag}.{out_cfg.fmt}"
            ),
        )

        # Figure 2b: m sweep at T/K = 2.0
        fig_m2 = fig_trajectories_m_sweep(
            data=data, task_maps=task_maps,
            data_cfg=data_cfg, fig_cfg=fig_cfg,
            target_tk=2.0,
            save_path=os.path.join(
                out_cfg.save_dir,
                f"{out_cfg.filename}_m_sweep_tk2.0{tag}.{out_cfg.fmt}"
            ),
        )

        if out_cfg.show:
            plt.show()
        else:
            plt.close(fig)
            for f_ in [fig_m1, fig_m2]:
                if f_ is not None:
                    plt.close(f_)