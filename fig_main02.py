#!/usr/bin/env python3
"""
fig_main.py
===========

Main figure for the paper.

Layout (1 × 2):
  Panel A: m = 1
  Panel B: m = T
  Both panels:
    x = T/K  (number of tasks scaled by number of programs)
    y = final realized phenotypic differentiation (normalized)
        i.e. mean across replicates of (final pheno_dist / task separation)
    One line-with-circles per dT value.
    Grayscale color ramp: low dT = light gray, high dT = black.
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
    densities: List[float] = field(default_factory=lambda: [0.1, 0.25, 0.5])
    T_values: List[int] = field(default_factory=lambda: [2, 4, 6, 8])
    task_divs: List[float] = field(default_factory=lambda: [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])


@dataclass
class OutputConfig:
    save_dir: str = "figures_out"
    filename: str = "fig_main"
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

def make_dt_gray_map(task_divs: List[float]) -> Dict[float, Tuple[float, float, float]]:
    """
    Map each dT to a grayscale RGB tuple.
    Uses linearly-spaced gray levels between 0.72 (lightest, still
    clearly visible on white) and 0.0 (black).  Equal perceptual
    spacing so every pair of adjacent lines is distinguishable.
    """
    vals = sorted(task_divs)
    if len(vals) == 1:
        return {vals[0]: (0.0, 0.0, 0.0)}

    gray_light = 0.72   # lightest line — visible but clearly light
    gray_dark = 0.0     # darkest line — black

    result = {}
    for i, dT in enumerate(vals):
        g = gray_light - (gray_light - gray_dark) * i / (len(vals) - 1)
        g = max(g, 0.0)   # clamp floating-point noise
        result[dT] = (g, g, g)
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
# 6. Final normalized pheno_dist helpers
# ============================================================

def compute_final_pheno_dist_norm(reps: List[Dict], target_dist: float) -> float:
    """
    Mean across replicates of the final normalized phenotypic distance.
    Each replicate's value = pheno_dist[final] / target_dist.
    """
    if not np.isfinite(target_dist) or target_dist <= 0:
        return np.nan
    vals = []
    for r in reps:
        n = r["n_actual_subs"]
        vals.append(float(r["pheno_dist"][n - 1]) / target_dist)
    return float(np.mean(vals)) if vals else np.nan


def compute_final_pheno_dist_norm_sem(reps: List[Dict], target_dist: float) -> float:
    """SEM across replicates."""
    if not np.isfinite(target_dist) or target_dist <= 0:
        return np.nan
    vals = []
    for r in reps:
        n = r["n_actual_subs"]
        vals.append(float(r["pheno_dist"][n - 1]) / target_dist)
    if len(vals) < 2:
        return 0.0
    return float(np.std(vals, ddof=1) / np.sqrt(len(vals)))


# ============================================================
# 6b. Angular projection error helpers
# ============================================================

def angular_projection_error(tasks: np.ndarray, phenotype: np.ndarray,
                             mode: str = "mean") -> float:
    """
    Angular mismatch between task targets and realized phenotypes.
    tasks: (L, T), phenotype: (L, T).
    Returns angle in radians.
    """
    vals = []
    for j in range(tasks.shape[1]):
        t = tasks[:, j]
        z = phenotype[:, j]
        nt = np.linalg.norm(t)
        nz = np.linalg.norm(z)
        if nt < 1e-12 or nz < 1e-12:
            vals.append(np.pi / 2)
            continue
        cos_th = np.dot(t, z) / (nt * nz)
        cos_th = np.clip(cos_th, -1.0, 1.0)
        vals.append(float(np.arccos(cos_th)))
    vals = np.array(vals)
    if mode == "mean":
        return float(np.mean(vals))
    if mode == "max":
        return float(np.max(vals))
    raise ValueError("mode must be 'mean' or 'max'")


def compute_final_angproj(reps: List[Dict], tasks: np.ndarray,
                          mode: str = "mean") -> float:
    """Mean angular projection error at final snapshot, averaged across replicates."""
    vals = []
    for r in reps:
        step = max(r["snapshots"].keys())
        phen = r["snapshots"][step]["phenotype"]
        vals.append(angular_projection_error(tasks, phen, mode=mode))
    return float(np.mean(vals)) if vals else np.nan


def compute_final_angproj_sem(reps: List[Dict], tasks: np.ndarray,
                              mode: str = "mean") -> float:
    """SEM of angular projection error across replicates."""
    vals = []
    for r in reps:
        step = max(r["snapshots"].keys())
        phen = r["snapshots"][step]["phenotype"]
        vals.append(angular_projection_error(tasks, phen, mode=mode))
    if len(vals) < 2:
        return 0.0
    return float(np.std(vals, ddof=1) / np.sqrt(len(vals)))


# ============================================================
# 7. MAIN FIGURE
# ============================================================

def fig_main(
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
      One line per dT, grayscale color ramp.
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

    dt_color_map = make_dt_gray_map(data_cfg.task_divs)

    # Precompute target distances for normalization
    target_dists = {
        T: {
            dT: mean_pairwise_distance(task_maps[T][dT])
            if (T in task_maps and dT in task_maps[T]) else np.nan
            for dT in data_cfg.task_divs
        }
        for T in data_cfg.T_values
    }

    # T/K x-values
    tk_values = np.array([T / data_cfg.K for T in data_cfg.T_values])

    panel_labels = "ABCDEFGH"

    for gi, group in enumerate(fig_cfg.groups):
        ax = axes[gi]
        ax.set_box_aspect(1)

        # Panel label
        ax.text(-0.15, 1.08, panel_labels[gi],
                transform=ax.transAxes, fontsize=14, fontweight="bold",
                va="top", ha="left")

        ax.set_title(group.label, fontsize=12, pad=6)

        # ── Line plot of final pheno_dist_norm vs T/K ──
        for dT in data_cfg.task_divs:
            color = dt_color_map[dT]
            y_means = []
            y_sems = []
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

                mean_val = compute_final_pheno_dist_norm(reps, td)
                sem_val = compute_final_pheno_dist_norm_sem(reps, td)

                if np.isfinite(mean_val):
                    x_plot.append(tk_values[ti])
                    y_means.append(mean_val)
                    y_sems.append(sem_val)

            if x_plot:
                x_plot = np.array(x_plot)
                y_means = np.array(y_means)
                y_sems = np.array(y_sems)

                ax.plot(x_plot, y_means, "-o", color=color,
                        lw=1.5, ms=5, label=f"ΔT = {dT}")
                ax.fill_between(x_plot, y_means - y_sems, y_means + y_sems,
                                color=color, alpha=0.12)

        # Format
        ax.axhline(1, color="gray", ls="--", lw=0.8, alpha=0.5)
        ax.set_ylim(0, 1.1)
        ax.set_xticks(tk_values)
        ax.set_xticklabels([f"{v:.1f}" for v in tk_values])
        ax.set_xlabel("T / K")
        if gi == 0:
            ax.set_ylabel("Realized differentiation\n(pheno dist / task separation)")
        else:
            ax.tick_params(labelleft=False)

    # Colorbar for dT
    _add_grayscale_colorbar(fig, data_cfg.task_divs)

    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
        print(f"Saved: {save_path}")

    return fig


# ============================================================
# 7b. ANGULAR PROJECTION ERROR FIGURE
# ============================================================

def fig_angproj(
    data: Dict,
    task_maps: Dict,
    data_cfg: DataConfig,
    fig_cfg: FigureConfig,
    mode: str = "mean",
    save_path: Optional[str] = None,
):
    """
    1 × 2 figure (same layout as fig_main):
      Panel A (m=1) and Panel B (m=T).
      x = T/K, y = mean angular projection error (radians).
      One line per dT, grayscale color ramp.
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

    dt_color_map = make_dt_gray_map(data_cfg.task_divs)
    tk_values = np.array([T / data_cfg.K for T in data_cfg.T_values])
    panel_labels = "ABCDEFGH"

    # Collect all y-values to set shared ylim
    all_y = []

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
            y_sems = []
            x_plot = []

            for ti, T in enumerate(data_cfg.T_values):
                if dT not in data[T] or not data[T][dT]:
                    continue
                if T not in task_maps or dT not in task_maps[T]:
                    continue

                reps_by_m = data[T][dT]
                available_m = sorted(reps_by_m.keys())
                m = resolve_m_selector(group.m_selector, available_m, T)
                if m is None or m not in reps_by_m:
                    continue

                reps = reps_by_m[m]
                tasks = task_maps[T][dT]

                mean_val = compute_final_angproj(reps, tasks, mode=mode)
                sem_val = compute_final_angproj_sem(reps, tasks, mode=mode)

                if np.isfinite(mean_val):
                    x_plot.append(tk_values[ti])
                    y_means.append(mean_val)
                    y_sems.append(sem_val)

            if x_plot:
                x_plot = np.array(x_plot)
                y_means = np.array(y_means)
                y_sems = np.array(y_sems)
                all_y.extend(y_means + y_sems)

                ax.plot(x_plot, y_means, "-o", color=color,
                        lw=1.5, ms=5, label=f"ΔT = {dT}")
                ax.fill_between(x_plot, y_means - y_sems, y_means + y_sems,
                                color=color, alpha=0.12)

        ax.set_xticks(tk_values)
        ax.set_xticklabels([f"{v:.1f}" for v in tk_values])
        ax.set_xlabel("T / K")
        if gi == 0:
            mode_label = "mean" if mode == "mean" else "max"
            ax.set_ylabel(f"Angular projection error\n({mode_label}, radians)")
        else:
            ax.tick_params(labelleft=False)

    # Shared y-axis: 0 to max + margin, identical across panels
    if all_y:
        y_top = max(all_y) * 1.1
    else:
        y_top = np.pi / 2
    for ax in axes:
        ax.set_ylim(0, y_top)
        ax.axhline(np.pi / 2, color="gray", ls="--", lw=0.8, alpha=0.5)

    # Colorbar for dT
    _add_grayscale_colorbar(fig, data_cfg.task_divs)

    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
        print(f"Saved: {save_path}")

    return fig


# ============================================================
# 8. COLORBAR
# ============================================================

def _add_grayscale_colorbar(fig, task_divs: List[float]):
    """Single grayscale colorbar for dT, centered below the panels."""
    dt_min = min(task_divs)
    dt_max = max(task_divs)
    dt_norm = mpl.colors.Normalize(vmin=dt_min, vmax=dt_max)

    # Light gray → black (matching make_dt_gray_map: 0.72 → 0.0)
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "dt_gray", [(0.72, 0.72, 0.72), (0.0, 0.0, 0.0)]
    )

    cb_w = 0.50
    cb_x = 0.5 - cb_w / 2 + 0.02
    cb_y = 0.02
    cb_h = 0.025

    ax_cb = fig.add_axes([cb_x, cb_y, cb_w, cb_h])
    cb = mpl.colorbar.ColorbarBase(
        ax_cb, cmap=cmap, norm=dt_norm,
        orientation="horizontal", ticks=task_divs
    )
    cb.set_ticklabels([str(dT) for dT in task_divs])
    cb.set_label(r"$\Delta T$", fontsize=10)
    cb.ax.tick_params(labelsize=8)


# ============================================================
# 9. NUMERICAL SUMMARY
# ============================================================

def print_numerical_summary(
    data: Dict,
    task_maps: Dict,
    data_cfg: DataConfig,
    fig_cfg: FigureConfig,
):
    print("\n" + "=" * 72)
    print("NUMERICAL SUMMARY")
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
        f"{'pheno_norm':>12}  {'sem':>10}"
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

                reps = reps_by_m[m]
                val = compute_final_pheno_dist_norm(reps, td)
                sem = compute_final_pheno_dist_norm_sem(reps, td)

                print(
                    f"{T / data_cfg.K:>6.1f}  {dT:>5.1f}  {group.label:>8}  {m:>4}  "
                    f"{val:>12.4f}  {sem:>10.4f}"
                )


def print_angproj_summary(
    data: Dict,
    task_maps: Dict,
    data_cfg: DataConfig,
    fig_cfg: FigureConfig,
    mode: str = "mean",
):
    mode_label = "mean" if mode == "mean" else "max"
    print(f"\n{'=' * 78}")
    print(f"ANGULAR PROJECTION ERROR SUMMARY  ({mode_label})")
    print("=" * 78)

    header = (
        f"{'T/K':>6}  {'dT':>5}  {'group':>8}  {'m':>4}  "
        f"{'ang_err':>12}  {'sem':>10}"
    )
    print(header)
    print("-" * len(header))

    for T in data_cfg.T_values:
        for dT in data_cfg.task_divs:
            if dT not in data[T] or not data[T][dT]:
                continue
            if T not in task_maps or dT not in task_maps[T]:
                continue

            reps_by_m = data[T][dT]
            available_m = sorted(reps_by_m.keys())
            tasks = task_maps[T][dT]

            for group in fig_cfg.groups:
                m = resolve_m_selector(group.m_selector, available_m, T)
                if m is None or m not in reps_by_m:
                    continue

                reps = reps_by_m[m]
                val = compute_final_angproj(reps, tasks, mode=mode)
                sem = compute_final_angproj_sem(reps, tasks, mode=mode)

                print(
                    f"{T / data_cfg.K:>6.1f}  {dT:>5.1f}  {group.label:>8}  {m:>4}  "
                    f"{val:>12.4f}  {sem:>10.4f}"
                )


# ============================================================
# 10. CLI
# ============================================================

def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument("--cache_dir", type=str, default="simulation_cache")
    p.add_argument("--save_dir", type=str, default="figures_out")
    p.add_argument("--filename", type=str, default="fig_main")
    p.add_argument("--fmt", type=str, default="pdf")

    p.add_argument("--L", type=int, default=100)
    p.add_argument("--gamma", type=float, default=1.0)
    p.add_argument(
        "--densities",
        type=float,
        nargs="+",
        default=[0.1, 0.25, 0.5],
        dest="densities")

    p.add_argument("--K", type=int, default=4)
    p.add_argument("--T", type=int, nargs="+", default=[2,4,6,8], dest="T_values")
    p.add_argument("--dT", type=float, nargs="+",
                   default=[0.2,0.4,0.6,0.8,1.0,1.2,1.4],
                   dest="task_divs")

    p.add_argument("--no_show", action="store_true")
    p.add_argument("--no_summary", action="store_true")

    return p.parse_args()


# ============================================================
# 11. MAIN
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

        print(f"Plotting density={density:.4f}...")
        fig = fig_main(
            data=data,
            task_maps=task_maps,
            data_cfg=data_cfg,
            fig_cfg=fig_cfg,
            save_path=save_path,
        )

        # Angular projection error figure
        angproj_path = os.path.join(
            out_cfg.save_dir,
            f"fig_angproj_density{density:.4f}.{out_cfg.fmt}"
        )
        print(f"Plotting angular projection error for density={density:.4f}...")
        fig_ang = fig_angproj(
            data=data,
            task_maps=task_maps,
            data_cfg=data_cfg,
            fig_cfg=fig_cfg,
            mode="mean",
            save_path=angproj_path,
        )

        if out_cfg.print_summary:
            print(f"\nSummary for density={density:.4f}")
            print_numerical_summary(data, task_maps, data_cfg, fig_cfg)
            print_angproj_summary(data, task_maps, data_cfg, fig_cfg, mode="mean")

        if out_cfg.show:
            plt.show()
        else:
            plt.close(fig)
            plt.close(fig_ang)