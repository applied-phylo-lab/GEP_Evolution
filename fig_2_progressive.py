#!/usr/bin/env python3
"""
fig_2_progressive.py
====================

Progressive-reveal version of fig_2 for presentations.

Plots ONLY the dT = 1.4 line (purple, darkest viridis_r color) for
pheno_dist_norm across T/K, with two panels (m=1 and m=T).

Generates a SEQUENCE of figures, one per cutoff in [1, 50, 100, 150, 200]:
  - Plot k shows cutoff[k] in full purple
  - All earlier cutoffs are drawn faintly in light gray (line + markers, low alpha)

Each figure is saved with a _step{k}_cutoff{c} tag so you can drop them
into slides in order for a build-up effect.

Only the main (pheno_dist) layout is produced. No usage figure, no colorbar
(only one dT is shown).
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
    densities: List[float] = field(default_factory=lambda: [1 / 4])
    T_values: List[int] = field(default_factory=lambda: [2, 3, 4, 6, 8])
    task_divs: List[float] = field(default_factory=lambda: [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])
    # The single dT we actually plot
    dT_focus: float = 1.4
    # Progressive cutoff sequence
    cutoffs: List[int] = field(default_factory=lambda: [1, 50, 100, 150, 200, 500])


@dataclass
class OutputConfig:
    save_dir: str = "figures_out"
    filename: str = "fig_2_progressive"
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
        GroupSpec(label="Sequential", m_selector="min"),
        GroupSpec(label="Simultaneous", m_selector="T"),
    ])
    # Styling
    faint_color: str = "lightgray"
    faint_alpha: float = 0.55
    # Y-axis limits held fixed across the progression so the sequence
    # of slides is visually comparable. Tuned for dT = 1.4 at the
    # default simulation parameters.
    ylim: Tuple[float, float] = (0.0, 1.15)


# ============================================================
# 2. HELPERS (subset of fig_2.py)
# ============================================================

def mean_pairwise_distance(X: np.ndarray) -> float:
    T = X.shape[1]
    if T < 2:
        return np.nan
    return float(np.mean([
        np.linalg.norm(X[:, i] - X[:, j])
        for i, j in combinations(range(T), 2)
    ]))


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


def pheno_dist_norm_at_step(reps: List[Dict], target_dist: float, cutoff: int) -> np.ndarray:
    """Normalized pheno_dist at cutoff. Stalled reps use their final value."""
    if not np.isfinite(target_dist) or target_dist <= 0:
        return np.array([])
    return np.array([
        float(r["pheno_dist"][min(cutoff - 1, r["n_actual_subs"] - 1)]) / target_dist
        for r in reps
    ])


def purple_color() -> Tuple:
    """The darkest viridis_r color — matches the dT=1.4 hue in the full figure."""
    return mpl.colormaps["viridis_r"](1.0)


# ============================================================
# 3. DATA LOADING (loads only dT_focus to save time)
# ============================================================

def load_focus_for_density(cfg: DataConfig, density: float):
    """Load only the single dT we care about, for each T."""
    data: Dict[int, Dict[float, Dict[int, List[Dict]]]] = {}
    task_maps: Dict[int, Dict[float, np.ndarray]] = {}

    for T in cfg.T_values:
        data[T] = {cfg.dT_focus: {}}
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

        if cfg.dT_focus not in alpha_map:
            print(f"WARNING: dT={cfg.dT_focus} not in alpha_map for T={T}")
            continue
        alpha = alpha_map[cfg.dT_focus]

        for m in range(1, T + 1):
            sp = sim_cache_path(cfg.cache_dir, cfg.L, cfg.K, cfg.gamma, cfg.fitness_r,
                                density, T, cfg.dT_focus, m, alpha)
            if os.path.exists(sp + ".npz"):
                results, _ = load_condition(sp)
                data[T][cfg.dT_focus][m] = results

        if data[T][cfg.dT_focus]:
            print(f"density={density:.4f}  T={T} dT={cfg.dT_focus}: "
                  f"m={sorted(data[T][cfg.dT_focus].keys())}")

    return data, task_maps


# ============================================================
# 4. PROGRESSIVE PLOT
# ============================================================

def _curve_for_cutoff(data, data_cfg, group, cutoff, target_dists):
    """
    Compute (x, mean, std) for the single dT_focus line at this cutoff.
    x-axis is T (number of tasks).
    """
    x_plot, y_means, y_stds = [], [], []

    dT = data_cfg.dT_focus
    for T in data_cfg.T_values:
        if T not in data or dT not in data[T] or not data[T][dT]:
            continue
        available_m = sorted(data[T][dT].keys())
        m = resolve_m_selector(group.m_selector, available_m, T)
        if m is None or m not in data[T][dT]:
            continue
        td = target_dists[T][dT]
        vals = pheno_dist_norm_at_step(data[T][dT][m], td, cutoff)
        mean_val, std_val = mean_and_std(vals)
        if np.isfinite(mean_val):
            x_plot.append(float(T))
            y_means.append(mean_val)
            y_stds.append(std_val)

    return np.array(x_plot), np.array(y_means), np.array(y_stds)


def fig_progressive(data, task_maps, data_cfg: DataConfig, fig_cfg: FigureConfig,
                    cutoff_idx: int, save_path: Optional[str] = None):
    """
    Build one figure of the progressive sequence.

    cutoff_idx is the index into data_cfg.cutoffs of the CURRENT (highlighted)
    cutoff. All cutoffs at indices < cutoff_idx are drawn faintly in gray.
    """
    current_cutoff = data_cfg.cutoffs[cutoff_idx]
    prior_cutoffs = data_cfg.cutoffs[:cutoff_idx]

    target_dists = {
        T: {dT: mean_pairwise_distance(task_maps[T][dT])
            if (T in task_maps and dT in task_maps[T]) else np.nan
            for dT in [data_cfg.dT_focus]}
        for T in data_cfg.T_values
    }

    n_groups = len(fig_cfg.groups)
    fig, axes = plt.subplots(1, n_groups, figsize=fig_cfg.figsize)
    if n_groups == 1:
        axes = np.array([axes])
    fig.subplots_adjust(wspace=0.25, left=0.12, right=0.92, top=0.92, bottom=0.15)

    t_values = np.array([float(T) for T in data_cfg.T_values])
    current_color = purple_color()

    for gi, group in enumerate(fig_cfg.groups):
        ax = axes[gi]
        ax.set_box_aspect(1)
        ax.set_title(group.label, fontsize=16, pad=8)

        # 1) Prior cutoffs as faint gray (drawn first, so they sit under current)
        for prior_c in prior_cutoffs:
            x, mu, sd = _curve_for_cutoff(data, data_cfg, group, prior_c, target_dists)
            if len(x) == 0:
                continue
            ax.errorbar(
                x, mu, yerr=sd,
                fmt="-o",
                color=fig_cfg.faint_color,
                alpha=fig_cfg.faint_alpha,
                lw=0.75, ms=5,
                markerfacecolor="none",
                markeredgecolor=fig_cfg.faint_color,
                capsize=2, capthick=1.0, elinewidth=1.0,
                zorder=2,
            )

        # 2) Current cutoff in full purple
        x, mu, sd = _curve_for_cutoff(data, data_cfg, group, current_cutoff, target_dists)
        if len(x) > 0:
            ax.errorbar(
                x, mu, yerr=sd,
                fmt="-o",
                color=current_color,
                lw=0.75, ms=5,
                markerfacecolor="none",
                markeredgecolor=current_color,
                capsize=2, capthick=1.0, elinewidth=1.0,
                zorder=3,
            )

        # Reference: dZ/dT = 1
        ax.axhline(1, color="gray", ls="--", lw=0.8, alpha=0.5)
        # Reference: T = K (number of tasks equals number of programs)
        ax.axvline(data_cfg.K, color="black", ls=":", lw=1.8, alpha=0.85, zorder=1)

        ax.set_ylim(*fig_cfg.ylim)
        ax.set_xticks(t_values)
        ax.set_xticklabels([f"{int(v)}" for v in t_values])
        ax.set_xlabel("Number of tasks", fontsize=14)
        if gi == 0:
            ax.set_ylabel(r"Phenotype divergence  $\bar{\Delta Z} / \bar{\Delta T}$",
                          fontsize=14)
        else:
            ax.tick_params(labelleft=False)

    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
        print(f"Saved: {save_path}")
    return fig


# ============================================================
# 5. NUMERICAL SUMMARY (optional, only for the dT_focus line)
# ============================================================

def print_summary(data, task_maps, data_cfg: DataConfig, fig_cfg: FigureConfig):
    target_dists = {
        T: {data_cfg.dT_focus: mean_pairwise_distance(task_maps[T][data_cfg.dT_focus])
            if (T in task_maps and data_cfg.dT_focus in task_maps[T]) else np.nan}
        for T in data_cfg.T_values
    }

    print(f"\n{'=' * 80}")
    print(f"SUMMARY — dT = {data_cfg.dT_focus}, progressive cutoffs")
    print("=" * 80)
    header = (f"{'cutoff':>7}  {'group':>8}  {'T/K':>5}  {'m':>3}  "
              f"{'pheno':>8}  {'std':>7}  {'n':>4}")
    print(header)
    print("-" * len(header))

    dT = data_cfg.dT_focus
    for cutoff in data_cfg.cutoffs:
        for group in fig_cfg.groups:
            for T in data_cfg.T_values:
                if T not in data or dT not in data[T] or not data[T][dT]:
                    continue
                available_m = sorted(data[T][dT].keys())
                m = resolve_m_selector(group.m_selector, available_m, T)
                if m is None or m not in data[T][dT]:
                    continue
                reps = data[T][dT][m]
                vals = pheno_dist_norm_at_step(reps, target_dists[T][dT], cutoff)
                mu, sd = mean_and_std(vals)
                print(f"{cutoff:>7d}  {group.label:>8}  "
                      f"{T / data_cfg.K:>5.1f}  {m:>3d}  "
                      f"{mu:>8.4f}  {sd:>7.4f}  {len(vals):>4d}")


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
    p.add_argument("--dT_focus",   type=float, default=d.dT_focus)
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
        dT_focus=args.dT_focus, cutoffs=args.cutoffs,
    )
    out_cfg = OutputConfig(
        save_dir=args.save_dir, filename=args.filename,
        fmt=args.fmt, show=not args.no_show,
        print_summary=not args.no_summary,
    )
    fig_cfg = FigureConfig()
    os.makedirs(out_cfg.save_dir, exist_ok=True)

    for density in data_cfg.densities:
        print(f"\nLoading data for density={density:.4f}  (dT={data_cfg.dT_focus} only)...")
        data, task_maps = load_focus_for_density(data_cfg, density)

        for k, cutoff in enumerate(data_cfg.cutoffs):
            sp = os.path.join(
                out_cfg.save_dir,
                f"{out_cfg.filename}_step{k}_sub{cutoff}"
                f"_dT{data_cfg.dT_focus}"
                f"_gamma{data_cfg.gamma}_fr{data_cfg.fitness_r}"
                f"_density{density:.4f}.{out_cfg.fmt}"
            )
            fig = fig_progressive(data, task_maps, data_cfg, fig_cfg,
                                  cutoff_idx=k, save_path=sp)
            if out_cfg.show:
                plt.show()
            else:
                plt.close(fig)

        if out_cfg.print_summary:
            print_summary(data, task_maps, data_cfg, fig_cfg)