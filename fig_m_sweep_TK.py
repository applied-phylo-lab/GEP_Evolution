#!/usr/bin/env python3
"""
fig_m_sweep_TK.py
=================

Same style as fig_main.py but sweeps all m values as columns.

Layout: 1 row x N_m columns (one per m value present in cache)
  x = T/K
  y = final normalized pheno_dist (mean pairwise dZ / dT)
  One line per dT, grayscale color ramp.
  ±std stick error bars.

For each column (m value), only T values where that m exists are plotted.
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
    densities: List[float] = field(default_factory=lambda: [0.25, 0.5])
    T_values: List[int] = field(default_factory=lambda: [2, 3, 4, 6, 8])
    task_divs: List[float] = field(default_factory=lambda: [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])


@dataclass
class OutputConfig:
    save_dir: str = "figures_out"
    filename: str = "fig_m_sweep_TK"
    fmt: str = "pdf"
    show: bool = True
    print_summary: bool = True


@dataclass
class FigureConfig:
    figsize_per_col: Tuple[float, float] = (4.0, 4.5)
    line_width: float = 0.75
    marker_size: float = 5.0


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


def compute_final_pheno_dist_norm_vals(reps: List[Dict], target_dist: float) -> np.ndarray:
    if not np.isfinite(target_dist) or target_dist <= 0:
        return np.array([])
    vals = []
    for r in reps:
        n = r["n_actual_subs"]
        vals.append(float(r["pheno_dist"][n - 1]) / target_dist)
    return np.array(vals)


def mean_and_std(vals: np.ndarray) -> Tuple[float, float]:
    if len(vals) == 0:
        return np.nan, np.nan
    if len(vals) == 1:
        return float(vals.mean()), 0.0
    return float(vals.mean()), float(vals.std(ddof=1))


# ============================================================
# 3. COLOR HELPERS
# ============================================================

def make_dt_gray_map(task_divs: List[float]) -> Dict[float, Tuple[float, float, float]]:
    vals = sorted(task_divs)
    if len(vals) == 1:
        return {vals[0]: (0.0, 0.0, 0.0)}
    gray_light = 0.72
    gray_dark = 0.0
    result = {}
    for i, dT in enumerate(vals):
        g = gray_light - (gray_light - gray_dark) * i / (len(vals) - 1)
        result[dT] = (max(g, 0.0),) * 3
    return result


# ============================================================
# 4. DATA LOADING
# ============================================================

def load_all_for_density(cfg: DataConfig, density: float):
    data = {}
    task_maps = {}
    alpha_maps = {}

    for T in cfg.T_values:
        data[T] = {}

        tpath = task_cache_path(cfg.cache_dir, cfg.L, cfg.K, cfg.gamma, T)
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
                    cfg.cache_dir, cfg.L, cfg.K, cfg.gamma,
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
# 5. FIGURE
# ============================================================

def fig_m_sweep_TK(
    data: Dict,
    task_maps: Dict,
    data_cfg: DataConfig,
    fig_cfg: FigureConfig,
    save_path: Optional[str] = None,
):
    """
    One column per m value (all m present across any condition).
    x = T/K, y = final pheno_dist_norm, lines = dT.
    """
    # Collect all m values present anywhere in the data
    all_m = sorted({
        m
        for T in data_cfg.T_values
        for dT in data_cfg.task_divs
        if T in data and dT in data[T]
        for m in data[T][dT].keys()
    })

    if not all_m:
        print("No data found.")
        return None

    n_cols = len(all_m)
    fig_w = fig_cfg.figsize_per_col[0] * n_cols
    fig_h = fig_cfg.figsize_per_col[1]

    fig, axes = plt.subplots(1, n_cols, figsize=(fig_w, fig_h))
    if n_cols == 1:
        axes = [axes]

    fig.subplots_adjust(
        wspace=0.15,
        left=0.10, right=0.96, top=0.88, bottom=0.15
    )

    dt_color_map = make_dt_gray_map(data_cfg.task_divs)
    tk_values = np.array([T / data_cfg.K for T in data_cfg.T_values])

    target_dists = {
        T: {
            dT: mean_pairwise_distance(task_maps[T][dT])
            if (T in task_maps and dT in task_maps[T]) else np.nan
            for dT in data_cfg.task_divs
        }
        for T in data_cfg.T_values
    }

    panel_labels = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    for col, m in enumerate(all_m):
        ax = axes[col]
        ax.set_box_aspect(1)

        ax.text(-0.15, 1.08, panel_labels[col],
                transform=ax.transAxes, fontsize=14, fontweight="bold",
                va="top", ha="left")
        ax.set_title(f"m = {m}", fontsize=12, pad=6)

        for dT in data_cfg.task_divs:
            color = dt_color_map[dT]
            x_plot = []
            y_means = []
            y_stds = []

            for ti, T in enumerate(data_cfg.T_values):
                if T not in data or dT not in data[T]:
                    continue
                if m not in data[T][dT]:
                    continue

                reps = data[T][dT][m]
                td = target_dists[T][dT]
                vals = compute_final_pheno_dist_norm_vals(reps, td)
                mean_val, std_val = mean_and_std(vals)

                if np.isfinite(mean_val):
                    x_plot.append(tk_values[ti])
                    y_means.append(mean_val)
                    y_stds.append(std_val)

            if x_plot:
                x_plot = np.array(x_plot)
                y_means = np.array(y_means)
                y_stds = np.array(y_stds)

                ax.errorbar(
                    x_plot, y_means, yerr=y_stds,
                    fmt="-o", color=color,
                    lw=fig_cfg.line_width, ms=fig_cfg.marker_size,
                    markerfacecolor="none", markeredgecolor=color,
                    capsize=2, capthick=1.0, elinewidth=1.0,
                )

        ax.axhline(1, color="gray", ls="--", lw=0.8, alpha=0.5)
        ax.set_ylim(0, None)
        ax.set_xticks(tk_values)
        ax.set_xticklabels([f"{v:.1f}" for v in tk_values])
        ax.set_xlabel("T / K")

        if col == 0:
            ax.set_ylabel(r"Mean pairwise $\bar{\Delta Z} / \bar{\Delta T}$")
        else:
            ax.tick_params(labelleft=False)

    _add_grayscale_colorbar(fig, data_cfg.task_divs)

    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
        print(f"Saved: {save_path}")

    return fig


# ============================================================
# 6. COLORBAR
# ============================================================

def _add_grayscale_colorbar(fig, task_divs: List[float]):
    dt_min = min(task_divs)
    dt_max = max(task_divs)
    dt_norm = mpl.colors.Normalize(vmin=dt_min, vmax=dt_max)
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "dt_gray", [(0.72, 0.72, 0.72), (0.0, 0.0, 0.0)]
    )
    cb_w = 0.50
    cb_x = 0.5 - cb_w / 2 + 0.02
    ax_cb = fig.add_axes([cb_x, 0.02, cb_w, 0.025])
    cb = mpl.colorbar.ColorbarBase(
        ax_cb, cmap=cmap, norm=dt_norm,
        orientation="horizontal", ticks=task_divs
    )
    cb.set_ticklabels([str(dT) for dT in task_divs])
    cb.set_label(r"$\Delta T$", fontsize=10)
    cb.ax.tick_params(labelsize=8)


# ============================================================
# 7. NUMERICAL SUMMARY
# ============================================================

def print_summary(data: Dict, task_maps: Dict, data_cfg: DataConfig):
    target_dists = {
        T: {
            dT: mean_pairwise_distance(task_maps[T][dT])
            if (T in task_maps and dT in task_maps[T]) else np.nan
            for dT in data_cfg.task_divs
        }
        for T in data_cfg.T_values
    }

    print("\n" + "=" * 72)
    print("NUMERICAL SUMMARY — final pheno_dist_norm")
    print("=" * 72)
    print(f"{'T/K':>6}  {'dT':>5}  {'m':>4}  {'mean':>10}  {'std':>10}")
    print("-" * 72)

    for T in data_cfg.T_values:
        for dT in data_cfg.task_divs:
            if T not in data or dT not in data[T] or not data[T][dT]:
                continue
            td = target_dists[T][dT]
            for m in sorted(data[T][dT].keys()):
                vals = compute_final_pheno_dist_norm_vals(data[T][dT][m], td)
                mean_val, std_val = mean_and_std(vals)
                print(f"{T/data_cfg.K:>6.1f}  {dT:>5.1f}  {m:>4}  "
                      f"{mean_val:>10.4f}  {std_val:>10.4f}")


# ============================================================
# 8. CLI
# ============================================================

def parse_args():
    defaults = DataConfig()
    p = argparse.ArgumentParser()
    p.add_argument("--cache_dir", type=str, default=defaults.cache_dir)
    p.add_argument("--save_dir", type=str, default="figures_out")
    p.add_argument("--filename", type=str, default="fig_m_sweep_TK")
    p.add_argument("--fmt", type=str, default="pdf")
    p.add_argument("--L", type=int, default=defaults.L)
    p.add_argument("--K", type=int, default=defaults.K)
    p.add_argument("--gamma", type=float, default=defaults.gamma)
    p.add_argument("--densities", type=float, nargs="+", default=defaults.densities)
    p.add_argument("--T", type=int, nargs="+", default=defaults.T_values, dest="T_values")
    p.add_argument("--dT", type=float, nargs="+", default=defaults.task_divs, dest="task_divs")
    p.add_argument("--no_show", action="store_true")
    p.add_argument("--no_summary", action="store_true")
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
        densities=args.densities,
        T_values=args.T_values,
        task_divs=args.task_divs,
    )

    out_cfg = OutputConfig(
        save_dir=args.save_dir,
        filename=args.filename,
        fmt=args.fmt,
        show=not args.no_show,
        print_summary=not args.no_summary,
    )

    fig_cfg = FigureConfig()
    os.makedirs(out_cfg.save_dir, exist_ok=True)

    for density in data_cfg.densities:
        print(f"\nLoading data for density={density:.4f}...")
        data, task_maps, _ = load_all_for_density(data_cfg, density)

        save_path = os.path.join(
            out_cfg.save_dir,
            f"{out_cfg.filename}_density{density:.4f}.{out_cfg.fmt}"
        )
        print(f"Plotting for density={density:.4f}...")
        fig = fig_m_sweep_TK(
            data=data,
            task_maps=task_maps,
            data_cfg=data_cfg,
            fig_cfg=fig_cfg,
            save_path=save_path,
        )

        if out_cfg.print_summary:
            print(f"\nSummary for density={density:.4f}")
            print_summary(data, task_maps, data_cfg)

        if out_cfg.show:
            plt.show()
        else:
            plt.close(fig)