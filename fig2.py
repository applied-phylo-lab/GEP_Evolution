#!/usr/bin/env python3
"""
fig2.py
=======

Main figure for the paper.

Layout (1 × 2):
  Panel A: m = 1
  Panel B: m = T
  Both panels:
    x = T/K  (number of tasks scaled by number of programs)
    y = final normalized pheno_dist and usage metrics
    One line-with-circles per dT value.
    Viridis color ramp: low dT = light, high dT = dark.
"""

import argparse
import json
import os
from dataclasses import dataclass, field
from itertools import combinations
from typing import Dict, List, Optional, Tuple, Union

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

MSelector = Union[int, str]   # int or one of {"min", "max", "T"}

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
    filename: str = "fig2"
    fmt: str = "pdf"
    show: bool = True
    print_summary: bool = True


@dataclass
class GroupSpec:
    label: str
    m_selector: MSelector


@dataclass
class FigureConfig:
    figsize: Tuple[float, float] = (9.0, 4.5)
    groups: List[GroupSpec] = field(default_factory=lambda: [
        GroupSpec(label="m = 1", m_selector="min"),
        GroupSpec(label="m = T", m_selector="T"),
    ])


# ============================================================
# 2. METRIC HELPERS
# ============================================================

def mean_pairwise_distance(X: np.ndarray) -> float:
    """Mean pairwise Euclidean distance between columns of X (L, T)."""
    T = X.shape[1]
    if T < 2:
        return np.nan
    return float(np.mean([
        np.linalg.norm(X[:, i] - X[:, j])
        for i, j in combinations(range(T), 2)
    ]))


# ============================================================
# 3. COLOR HELPERS
# ============================================================

def make_dt_viridis_map(task_divs: List[float]) -> Dict[float, Tuple]:
    """Map each dT to a viridis color, linearly spaced over the dT range."""
    vals = sorted(task_divs)
    cmap = mpl.colormaps["viridis_r"]
    if len(vals) == 1:
        return {vals[0]: cmap(0.5)}
    result = {}
    for i, dT in enumerate(vals):
        result[dT] = cmap(i / (len(vals) - 1))
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
# 5. GROUP / CONDITION HELPERS
# ============================================================

def resolve_m_selector(m_selector: MSelector, available_m: List[int], T: int) -> Optional[int]:
    if not available_m:
        return None
    if isinstance(m_selector, int):
        return m_selector if m_selector in available_m else None
    if m_selector == "min":
        return min(available_m)
    if m_selector == "max":
        return max(available_m)
    if m_selector == "T":
        return T if T in available_m else None
    raise ValueError(f"Unknown m_selector: {m_selector}")


# ============================================================
# 6. METRIC COMPUTATION
# ============================================================

def compute_final_pheno_dist_norm_vals(reps: List[Dict], target_dist: float) -> np.ndarray:
    """Per-replicate final normalized phenotypic distance."""
    if not np.isfinite(target_dist) or target_dist <= 0:
        return np.array([])
    vals = []
    for r in reps:
        n = r["n_actual_subs"]
        vals.append(float(r["pheno_dist"][n - 1]) / target_dist)
    return np.array(vals)


def cosine_sim_usage(usage: np.ndarray) -> float:
    """Mean pairwise cosine similarity across task usage vectors. usage: (T, K)."""
    A = np.array(usage, dtype=float)
    norms = np.linalg.norm(A, axis=1, keepdims=True)
    norms = np.where(norms < 1e-12, 1.0, norms)
    A = A / norms
    vals = []
    for i, j in combinations(range(A.shape[0]), 2):
        vals.append(float(np.dot(A[i], A[j])))
    return float(np.mean(vals)) if vals else 0.0


def mean_usage_concentration(usage: np.ndarray) -> float:
    """Mean per-task program concentration (1 - normalized entropy). usage: (T, K).
    1 = task uses exactly one program, 0 = task uses all programs equally."""
    A = np.array(usage, dtype=float)
    K = A.shape[1]
    row_sums = A.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums < 1e-12, 1.0, row_sums)
    P = A / row_sums
    with np.errstate(divide="ignore", invalid="ignore"):
        logP = np.where(P > 0, np.log(P), 0.0)
    H = -np.sum(P * logP, axis=1)
    if K > 1:
        H = H / np.log(K)
    return float(np.mean(1.0 - H))


def compute_final_usage_metric_vals(
    reps: List[Dict],
    metric: str = "cosine",
) -> np.ndarray:
    """Per-replicate final usage metric."""
    vals = []
    for r in reps:
        step = max(r["snapshots"].keys())
        usage = r["snapshots"][step]["usage"]
        if metric == "cosine":
            vals.append(cosine_sim_usage(usage))
        elif metric == "concentration":
            vals.append(mean_usage_concentration(usage))
        else:
            raise ValueError("metric must be 'cosine' or 'concentration'")
    return np.array(vals)


def mean_and_std(vals: np.ndarray) -> Tuple[float, float]:
    """Return mean and std with safe behavior for small n."""
    if len(vals) == 0:
        return np.nan, np.nan
    if len(vals) == 1:
        return float(vals.mean()), 0.0
    return float(vals.mean()), float(vals.std(ddof=1))


# ============================================================
# 7. MAIN FIGURE
# ============================================================

def fig_main_plot(
    data: Dict,
    task_maps: Dict,
    data_cfg: DataConfig,
    fig_cfg: FigureConfig,
    save_path: Optional[str] = None,
):
    """
    1 × 2 figure:
      Panel A (m=1) and Panel B (m=T).
      x = T/K, y = final normalized pheno_dist.
      One line per dT, grayscale color ramp, ±std error bars.
    """
    n_groups = len(fig_cfg.groups)
    fig_w, fig_h = fig_cfg.figsize
    fig, axes = plt.subplots(1, n_groups, figsize=(fig_w, fig_h))
    if n_groups == 1:
        axes = np.array([axes])
    fig.subplots_adjust(
        wspace=0.25,
        left=0.12, right=0.92, top=0.88, bottom=0.15
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

    tk_values = np.array([T / data_cfg.K for T in data_cfg.T_values])
    panel_labels = "ABCDEFGH"

    for gi, group in enumerate(fig_cfg.groups):
        ax = axes[gi]
        ax.set_box_aspect(1)

        ax.text(-0.15, 1.08, panel_labels[gi],
                transform=ax.transAxes, fontsize=14, fontweight="bold",
                va="top", ha="left")
        ax.set_title(group.label, fontsize=12, pad=6)

        for dT in data_cfg.task_divs:
            color = dt_color_map[dT]
            y_means = []
            y_stds = []
            x_plot = []

            for ti, T in enumerate(data_cfg.T_values):
                if dT not in data[T] or not data[T][dT]:
                    continue

                reps_by_m = data[T][dT]
                available_m = sorted(reps_by_m.keys())
                m = resolve_m_selector(group.m_selector, available_m, T)
                if m is None or m not in reps_by_m:
                    continue

                reps = reps_by_m[m]
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
                    lw=0.75, ms=5,
                    markerfacecolor="none", markeredgecolor=color,
                    capsize=2, capthick=1.0, elinewidth=1.0,
                )

        ax.axhline(1, color="gray", ls="--", lw=0.8, alpha=0.5)
        ax.set_ylim(0, None)
        ax.set_xticks(tk_values)
        ax.set_xticklabels([f"{v:.1f}" for v in tk_values])
        ax.set_xlabel("T / K")
        if gi == 0:
            ax.set_ylabel(r"Mean pairwise $\bar{\Delta Z} / \bar{\Delta T}$")
        else:
            ax.tick_params(labelleft=False)

    _add_viridis_colorbar(fig, data_cfg.task_divs)

    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
        print(f"Saved: {save_path}")

    return fig


# ============================================================
# 8. USAGE FIGURE
# ============================================================

def fig_usage_subset(
    data: Dict,
    data_cfg: DataConfig,
    fig_cfg: FigureConfig,
    save_path: Optional[str] = None,
):
    """
    2 × 2 figure:
      columns = m=1, m=T
      rows    = task usage similarity, program usage concentration
      x = T/K, lines = dT values — both driven entirely by data_cfg.
      ±std error bars.
    """
    regimes = fig_cfg.groups
    if len(regimes) != 2:
        raise ValueError("fig_usage_subset assumes exactly two groups: m=1 and m=T")

    fig, axes = plt.subplots(2, 2, figsize=(9.0, 8.0))
    fig.subplots_adjust(
        hspace=0.30, wspace=0.22,
        left=0.12, right=0.94, top=0.90, bottom=0.12
    )

    dt_color_map = make_dt_viridis_map(data_cfg.task_divs)
    tk_values = np.array([T / data_cfg.K for T in data_cfg.T_values], dtype=float)
    panel_labels = "ABCDEFGH"

    row_specs = [
        ("cosine",       r"Mean pairwise program usage similarity"),
        ("concentration", r"Mean per-task program usage concentration"),
    ]

    for col, group in enumerate(regimes):
        for row, (metric, ylabel) in enumerate(row_specs):
            ax = axes[row, col]
            ax.set_box_aspect(1)

            ax.text(-0.15, 1.08, panel_labels[row * 2 + col],
                    transform=ax.transAxes, fontsize=14, fontweight="bold",
                    va="top", ha="left")

            if row == 0:
                ax.set_title(group.label, fontsize=12, pad=6)

            for dT in data_cfg.task_divs:
                color = dt_color_map[dT]
                x_plot = []
                y_means = []
                y_stds = []

                for T in data_cfg.T_values:
                    if T not in data or dT not in data[T] or not data[T][dT]:
                        continue

                    reps_by_m = data[T][dT]
                    available_m = sorted(reps_by_m.keys())
                    m = resolve_m_selector(group.m_selector, available_m, T)
                    if m is None or m not in reps_by_m:
                        continue

                    reps = reps_by_m[m]
                    vals = compute_final_usage_metric_vals(reps, metric=metric)
                    mean_val, std_val = mean_and_std(vals)

                    if np.isfinite(mean_val):
                        x_plot.append(T / data_cfg.K)
                        y_means.append(mean_val)
                        y_stds.append(std_val)

                if x_plot:
                    x_plot = np.array(x_plot)
                    y_means = np.array(y_means)
                    y_stds = np.array(y_stds)

                    ax.errorbar(
                        x_plot, y_means, yerr=y_stds,
                        fmt="-o", color=color,
                        lw=0.75, ms=5,
                        markerfacecolor="none", markeredgecolor=color,
                        capsize=2, capthick=1.0, elinewidth=1.0,
                    )

            ax.set_xticks(tk_values)
            ax.set_xticklabels([f"{x:.1f}" for x in tk_values])
            ax.set_xlabel("T / K")

            if col == 0:
                ax.set_ylabel(ylabel)
            else:
                ax.tick_params(labelleft=False)

            ax.set_ylim(0, 1.05)

    _add_viridis_colorbar(fig, data_cfg.task_divs)

    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
        print(f"Saved: {save_path}")

    return fig


# ============================================================
# 9. COLORBAR
# ============================================================

def _add_viridis_colorbar(fig, task_divs: List[float]):
    """Single viridis colorbar for dT, centered below the panels."""
    dt_norm = mpl.colors.Normalize(vmin=min(task_divs), vmax=max(task_divs))
    cmap = mpl.colormaps["viridis_r"]
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
# 10. NUMERICAL SUMMARIES
# ============================================================

def print_numerical_summary(
    data: Dict,
    task_maps: Dict,
    data_cfg: DataConfig,
    fig_cfg: FigureConfig,
):
    print("\n" + "=" * 72)
    print("NUMERICAL SUMMARY — pheno_dist_norm")
    print("=" * 72)

    target_dists = {
        T: {
            dT: mean_pairwise_distance(task_maps[T][dT])
            if (T in task_maps and dT in task_maps[T]) else np.nan
            for dT in data_cfg.task_divs
        }
        for T in data_cfg.T_values
    }

    header = (
        f"{'T/K':>6}  {'dT':>5}  {'group':>8}  {'m':>4}  "
        f"{'mean':>10}  {'std':>10}"
    )
    print(header)
    print("-" * len(header))

    for T in data_cfg.T_values:
        for dT in data_cfg.task_divs:
            if dT not in data[T] or not data[T][dT]:
                continue

            reps_by_m = data[T][dT]
            available_m = sorted(reps_by_m.keys())
            td = target_dists[T][dT]

            for group in fig_cfg.groups:
                m = resolve_m_selector(group.m_selector, available_m, T)
                if m is None or m not in reps_by_m:
                    continue

                vals = compute_final_pheno_dist_norm_vals(reps_by_m[m], td)
                mean_val, std_val = mean_and_std(vals)

                print(
                    f"{T / data_cfg.K:>6.1f}  {dT:>5.1f}  {group.label:>8}  {m:>4}  "
                    f"{mean_val:>10.4f}  {std_val:>10.4f}"
                )


def print_usage_summary(
    data: Dict,
    data_cfg: DataConfig,
    fig_cfg: FigureConfig,
):
    print(f"\n{'=' * 78}")
    print("USAGE SUMMARY")
    print("=" * 78)

    header = (
        f"{'T/K':>6}  {'dT':>5}  {'group':>8}  {'m':>4}  "
        f"{'cosine_sim':>10}  {'cos_std':>10}  "
        f"{'concentr':>10}  {'con_std':>10}"
    )
    print(header)
    print("-" * len(header))

    for T in data_cfg.T_values:
        for dT in data_cfg.task_divs:
            if T not in data or dT not in data[T] or not data[T][dT]:
                continue

            reps_by_m = data[T][dT]
            available_m = sorted(reps_by_m.keys())

            for group in fig_cfg.groups:
                m = resolve_m_selector(group.m_selector, available_m, T)
                if m is None or m not in reps_by_m:
                    continue

                reps = reps_by_m[m]

                cos_vals = compute_final_usage_metric_vals(reps, metric="cosine")
                ent_vals = compute_final_usage_metric_vals(reps, metric="concentration")
                cos_mean, cos_std = mean_and_std(cos_vals)
                ent_mean, ent_std = mean_and_std(ent_vals)

                print(
                    f"{T / data_cfg.K:>6.1f}  {dT:>5.1f}  {group.label:>8}  {m:>4}  "
                    f"{cos_mean:>10.4f}  {cos_std:>10.4f}  "
                    f"{ent_mean:>10.4f}  {ent_std:>10.4f}"
                )


# ============================================================
# 11. CLI
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
    p.add_argument("--densities",  type=float, nargs="+", default=d_cfg.densities)
    p.add_argument("--T",          type=int,   nargs="+", default=d_cfg.T_values,  dest="T_values")
    p.add_argument("--dT",         type=float, nargs="+", default=d_cfg.task_divs, dest="task_divs")
    p.add_argument("--no_show",    action="store_true")
    p.add_argument("--no_summary", action="store_true")

    return p.parse_args()


# ============================================================
# 12. MAIN
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

        # Figure 1: main outcome
        save_path = os.path.join(
            out_cfg.save_dir,
            f"{out_cfg.filename}_density{density:.4f}.{out_cfg.fmt}"
        )
        print(f"Plotting main figure for density={density:.4f}...")
        fig = fig_main_plot(
            data=data,
            task_maps=task_maps,
            data_cfg=data_cfg,
            fig_cfg=fig_cfg,
            save_path=save_path,
        )

        # Figure 2: usage mechanism
        usage_path = os.path.join(
            out_cfg.save_dir,
            f"{out_cfg.filename}_usage_density{density:.4f}.{out_cfg.fmt}"
        )
        print(f"Plotting usage figure for density={density:.4f}...")
        fig_usage = fig_usage_subset(
            data=data,
            data_cfg=data_cfg,
            fig_cfg=fig_cfg,
            save_path=usage_path,
        )

        if out_cfg.print_summary:
            print(f"\nSummary for density={density:.4f}")
            print_numerical_summary(data, task_maps, data_cfg, fig_cfg)
            print_usage_summary(data, data_cfg, fig_cfg)

        if out_cfg.show:
            plt.show()
        else:
            plt.close(fig)
            plt.close(fig_usage)