#!/usr/bin/env python3
"""
fig_2.py
========

Same layout as fig2.py but plots pheno_dist and usage metrics at
arbitrary early substitution cutoffs rather than the final state.
Replicates that stalled before the cutoff use their actual final value.
Usage metrics use the snapshot closest to (but not exceeding) the cutoff.

Produces two figures per (density, cutoff): main and usage.
Pass --cutoffs 50 100 150 200 to inspect any set of cutoffs.
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

MSelector = Union[int, str]


@dataclass
class DataConfig:
    cache_dir: str = "simulation_cache"
    L: int = 100
    K: int = 4
    gamma: float = 1.0
    fitness_r: float = 0.0
    densities: List[float] = field(default_factory=lambda: [1/4])
    T_values: List[int] = field(default_factory=lambda: [2, 3, 4, 6, 8])
    task_divs: List[float] = field(default_factory=lambda: [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])
    cutoffs: List[int] = field(default_factory=lambda: [200])



@dataclass
class OutputConfig:
    save_dir: str = "figures_out"
    filename: str = "fig_2"
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


def make_dt_viridis_map(task_divs: List[float]) -> Dict[float, Tuple]:
    vals = sorted(task_divs)
    cmap = mpl.colormaps["viridis_r"]
    if len(vals) == 1:
        return {vals[0]: cmap(0.5)}
    return {dT: cmap(i / (len(vals) - 1)) for i, dT in enumerate(vals)}


def _add_viridis_colorbar(fig, task_divs: List[float]):
    dt_norm = mpl.colors.Normalize(vmin=min(task_divs), vmax=max(task_divs))
    cmap = mpl.colormaps["viridis_r"]
    cb_w = 0.50
    ax_cb = fig.add_axes([0.5 - cb_w / 2 + 0.02, 0.02, cb_w, 0.025])
    cb = mpl.colorbar.ColorbarBase(
        ax_cb, cmap=cmap, norm=dt_norm,
        orientation="horizontal", ticks=task_divs
    )
    cb.set_ticklabels([str(dT) for dT in task_divs])
    cb.set_label(r"$\Delta T$", fontsize=10)
    cb.ax.tick_params(labelsize=8)


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


def mean_and_std(vals: np.ndarray) -> Tuple[float, float]:
    if len(vals) == 0:
        return np.nan, np.nan
    if len(vals) == 1:
        return float(vals.mean()), 0.0
    return float(vals.mean()), float(vals.std(ddof=1))


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
                continue
            alpha = alpha_map[dT]
            for m in range(1, T + 1):
                sp = sim_cache_path(cfg.cache_dir, cfg.L, cfg.K, cfg.gamma, cfg.fitness_r,
                                    density, T, dT, m, alpha)
                if os.path.exists(sp + ".npz"):
                    results, _ = load_condition(sp)
                    data[T][dT][m] = results

        for dT in cfg.task_divs:
            if data[T].get(dT):
                print(f"density={density:.4f}  T={T} dT={dT}: m={sorted(data[T][dT].keys())}")

    return data, task_maps, alpha_maps


# ============================================================
# 4. METRIC COMPUTATION
# ============================================================

def pheno_dist_norm_at_step(reps: List[Dict], target_dist: float, cutoff: int) -> np.ndarray:
    """Normalized pheno_dist at cutoff. Stalled reps use their final value."""
    if not np.isfinite(target_dist) or target_dist <= 0:
        return np.array([])
    return np.array([
        float(r["pheno_dist"][min(cutoff - 1, r["n_actual_subs"] - 1)]) / target_dist
        for r in reps
    ])


def cosine_sim_usage(usage: np.ndarray) -> float:
    A = np.array(usage, dtype=float)
    norms = np.linalg.norm(A, axis=1, keepdims=True)
    A = A / np.where(norms < 1e-12, 1.0, norms)
    pairs = list(combinations(range(A.shape[0]), 2))
    return float(np.mean([np.dot(A[i], A[j]) for i, j in pairs])) if pairs else 0.0


def mean_usage_concentration(usage: np.ndarray) -> float:
    A = np.array(usage, dtype=float)
    K = A.shape[1]
    row_sums = A.sum(axis=1, keepdims=True)
    P = A / np.where(row_sums < 1e-12, 1.0, row_sums)
    with np.errstate(divide="ignore", invalid="ignore"):
        logP = np.where(P > 0, np.log(P), 0.0)
    H = -np.sum(P * logP, axis=1)
    if K > 1:
        H = H / np.log(K)
    return float(np.mean(1.0 - H))


def usage_metric_at_step(reps: List[Dict], cutoff: int, metric: str) -> np.ndarray:
    """Usage metric at snapshot closest to but not exceeding cutoff."""
    fn = cosine_sim_usage if metric == "cosine" else mean_usage_concentration
    vals = []
    for r in reps:
        snap_steps = sorted(r["snapshots"].keys())
        eligible = [s for s in snap_steps if s <= cutoff - 1]
        step = eligible[-1] if eligible else snap_steps[0]
        vals.append(fn(r["snapshots"][step]["usage"]))
    return np.array(vals)


# ============================================================
# 5. SHARED PLOT LOOP
# ============================================================

def _plot_lines(ax, data, data_cfg, fig_cfg, group, dt_color_map, get_vals_fn):
    """Plot errorbar lines over T/K for all dT into ax."""
    tk_values = np.array([T / data_cfg.K for T in data_cfg.T_values])

    for dT in data_cfg.task_divs:
        color = dt_color_map[dT]
        x_plot, y_means, y_stds = [], [], []

        for ti, T in enumerate(data_cfg.T_values):
            if dT not in data[T] or not data[T][dT]:
                continue
            m = resolve_m_selector(group.m_selector, sorted(data[T][dT].keys()), T)
            if m is None or m not in data[T][dT]:
                continue
            mean_val, std_val = mean_and_std(get_vals_fn(data[T][dT][m], T, dT))
            if np.isfinite(mean_val):
                x_plot.append(tk_values[ti])
                y_means.append(mean_val)
                y_stds.append(std_val)

        if x_plot:
            ax.errorbar(
                np.array(x_plot), np.array(y_means), yerr=np.array(y_stds),
                fmt="-o", color=color, lw=0.75, ms=5,
                markerfacecolor="none", markeredgecolor=color,
                capsize=2, capthick=1.0, elinewidth=1.0,
            )

    ax.set_xticks(tk_values)
    ax.set_xticklabels([f"{v:.1f}" for v in tk_values])
    ax.set_xlabel("T / K")


# ============================================================
# 6. FIGURES
# ============================================================

def fig_main_plot(data, task_maps, data_cfg, fig_cfg, cutoff, save_path=None):
    """1 × 2: pheno_dist_norm at cutoff."""
    target_dists = {
        T: {dT: mean_pairwise_distance(task_maps[T][dT])
            if (T in task_maps and dT in task_maps[T]) else np.nan
            for dT in data_cfg.task_divs}
        for T in data_cfg.T_values
    }

    n_groups = len(fig_cfg.groups)
    fig, axes = plt.subplots(1, n_groups, figsize=fig_cfg.figsize)
    if n_groups == 1:
        axes = np.array([axes])
    fig.subplots_adjust(wspace=0.25, left=0.12, right=0.92, top=0.88, bottom=0.15)

    dt_color_map = make_dt_viridis_map(data_cfg.task_divs)
    panel_labels = "ABCDEFGH"

    for gi, group in enumerate(fig_cfg.groups):
        ax = axes[gi]
        ax.set_box_aspect(1)
        ax.text(-0.15, 1.08, panel_labels[gi], transform=ax.transAxes,
                fontsize=14, fontweight="bold", va="top", ha="left")
        ax.set_title(f"{group.label}  (sub $\\leq$ {cutoff})", fontsize=12, pad=6)

        _plot_lines(ax, data, data_cfg, fig_cfg, group, dt_color_map,
                    lambda reps, T, dT: pheno_dist_norm_at_step(reps, target_dists[T][dT], cutoff))

        ax.axhline(1, color="gray", ls="--", lw=0.8, alpha=0.5)
        ax.set_ylim(0, None)
        if gi == 0:
            ax.set_ylabel(r"Mean pairwise $\bar{\Delta Z} / \bar{\Delta T}$")
        else:
            ax.tick_params(labelleft=False)

    _add_viridis_colorbar(fig, data_cfg.task_divs)
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
        print(f"Saved: {save_path}")
    return fig


def fig_usage_subset(data, data_cfg, fig_cfg, cutoff, save_path=None):
    """2 × 2: usage metrics at snapshot closest to cutoff."""
    if len(fig_cfg.groups) != 2:
        raise ValueError("fig_usage_subset assumes exactly two groups")

    row_specs = [
        ("cosine",        r"Mean pairwise program usage similarity"),
        ("concentration", r"Mean per-task program usage concentration"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(9.0, 8.0))
    fig.subplots_adjust(hspace=0.30, wspace=0.22,
                        left=0.12, right=0.94, top=0.90, bottom=0.12)

    dt_color_map = make_dt_viridis_map(data_cfg.task_divs)
    panel_labels = "ABCD"

    for col, group in enumerate(fig_cfg.groups):
        for row, (metric, ylabel) in enumerate(row_specs):
            ax = axes[row, col]
            ax.set_box_aspect(1)
            ax.text(-0.15, 1.08, panel_labels[row * 2 + col],
                    transform=ax.transAxes, fontsize=14, fontweight="bold",
                    va="top", ha="left")
            if row == 0:
                ax.set_title(f"{group.label}  (sub $\\leq$ {cutoff})", fontsize=12, pad=6)

            _plot_lines(ax, data, data_cfg, fig_cfg, group, dt_color_map,
                        lambda reps, T, dT, m=metric: usage_metric_at_step(reps, cutoff, m))

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
# 7. NUMERICAL SUMMARY
# ============================================================

def print_summary(data, task_maps, data_cfg, fig_cfg, cutoff):
    target_dists = {
        T: {dT: mean_pairwise_distance(task_maps[T][dT])
            if (T in task_maps and dT in task_maps[T]) else np.nan
            for dT in data_cfg.task_divs}
        for T in data_cfg.T_values
    }

    print(f"\n{'=' * 88}")
    print(f"SUMMARY — sub ≤ {cutoff}")
    print("=" * 88)
    header = (f"{'T/K':>6}  {'dT':>5}  {'group':>8}  {'m':>4}  "
              f"{'pheno':>10}  {'ph_std':>8}  "
              f"{'cos_sim':>8}  {'co_std':>7}  "
              f"{'conc':>8}  {'cn_std':>7}")
    print(header)
    print("-" * len(header))

    for T in data_cfg.T_values:
        for dT in data_cfg.task_divs:
            if dT not in data[T] or not data[T][dT]:
                continue
            td = target_dists[T][dT]
            available_m = sorted(data[T][dT].keys())

            for group in fig_cfg.groups:
                m = resolve_m_selector(group.m_selector, available_m, T)
                if m is None or m not in data[T][dT]:
                    continue
                reps = data[T][dT][m]

                ph_m, ph_s = mean_and_std(pheno_dist_norm_at_step(reps, td, cutoff))
                co_m, co_s = mean_and_std(usage_metric_at_step(reps, cutoff, "cosine"))
                cn_m, cn_s = mean_and_std(usage_metric_at_step(reps, cutoff, "concentration"))

                print(f"{T / data_cfg.K:>6.1f}  {dT:>5.1f}  {group.label:>8}  {m:>4}  "
                      f"{ph_m:>10.4f}  {ph_s:>8.4f}  "
                      f"{co_m:>8.4f}  {co_s:>7.4f}  "
                      f"{cn_m:>8.4f}  {cn_s:>7.4f}")


# ============================================================
# 8. CLI + MAIN
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
    fig_cfg = FigureConfig()
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
                    f"_density{density:.4f}.{out_cfg.fmt}"
                )

            fig = fig_main_plot(data, task_maps, data_cfg, fig_cfg,
                                cutoff=cutoff, save_path=sp("main"))
            fig_u = fig_usage_subset(data, data_cfg, fig_cfg,
                                     cutoff=cutoff, save_path=sp("usage"))

            if out_cfg.print_summary:
                print_summary(data, task_maps, data_cfg, fig_cfg, cutoff)

            if out_cfg.show:
                plt.show()
            else:
                plt.close(fig)
                plt.close(fig_u)