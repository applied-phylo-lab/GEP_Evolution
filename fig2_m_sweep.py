#!/usr/bin/env python3
"""
fig2_m_sweep.py
===============

Flexible plotting for cached simulation outputs.

Default behavior:
- One figure per density
- 1 row x N columns, one panel per T/K ratio
- x-axis: m (simultaneity, 1 to T)
- y-axis: endpoint pheno_dist / mean pairwise task-target distance
- lines: one per dT
- shaded band: mean ± SEM across replicates
- optional numerical summary printing

Key improvements:
- parameters are centralized in dataclasses
- density can be a list
- panel count is determined automatically from T_values
- dT colors are generalized instead of hard-coded to only 3 values
- output filenames automatically include density
- missing/invalid conditions are handled more explicitly
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
    densities: List[float] = field(default_factory=lambda: [0.1, 0.25, 0.5])
    T_values: List[int] = field(default_factory=lambda: [2, 4, 6, 8])
    task_divs: List[float] = field(default_factory=lambda: [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])

@dataclass
class OutputConfig:
    save_dir: str = "figures_out"
    filename: str = "fig2_m_sweep"
    fmt: str = "pdf"
    show: bool = True
    print_summary: bool = True


@dataclass
class FigureConfig:
    figsize_per_col: Tuple[float, float] = (4.5, 4.5)
    line_width: float = 1.8
    marker_size: float = 5.0
    band_alpha: float = 0.15
    reference_line_y: float = 1.0


DEFAULT_DT_COLORS = {
    0.2: "#4878CF",
    0.6: "#6ACC65",
    1.2: "#D65F5F",
}


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


def endpoint_pheno_dist_norm(reps: List[Dict], target_dist: float) -> np.ndarray:
    """
    Endpoint normalized phenotype distance for each replicate.
    Returns an array of per-replicate endpoint values.
    """
    if not (np.isfinite(target_dist) and target_dist > 0):
        return np.array([], dtype=float)

    vals = []
    for r in reps:
        n = r["n_actual_subs"]
        if n <= 0:
            continue
        vals.append(float(r["pheno_dist"][n - 1]) / target_dist)
    return np.array(vals, dtype=float)


def mean_and_sem(vals: np.ndarray) -> Tuple[float, float]:
    """Return mean and SEM with safe behavior for small n."""
    if len(vals) == 0:
        return np.nan, np.nan
    if len(vals) == 1:
        return float(vals.mean()), 0.0
    return float(vals.mean()), float(vals.std(ddof=1) / np.sqrt(len(vals)))


# ============================================================
# 3. COLOR HELPERS
# ============================================================

def make_dt_color_map(task_divs: List[float]) -> Dict[float, Tuple[float, float, float, float]]:
    """
    Map dT values to grayscale colors from white (low) to black (high),
    normalized over [0, sqrt(2)].
    """
    dT_min, dT_max = 0.0, np.sqrt(2)

    def normalize(dT):
        return (dT - dT_min) / (dT_max - dT_min)

    color_map = {}
    for dT in task_divs:
        x = normalize(dT)
        x = np.clip(x, 0.0, 1.0)

        # grayscale: (1-x) gives white→black
        gray = 1.0 - x
        color_map[dT] = (gray, gray, gray, 1.0)

    return color_map


def make_dt_label(dT: float) -> str:
    return rf"$\Delta T = {dT}$"


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
                print(f"WARNING: dT={dT} missing from alpha_map for T={T}")
                continue

            alpha = alpha_map[dT]
            missing_ms = []

            for m in range(1, T + 1):
                sp = sim_cache_path(
                    cfg.cache_dir, cfg.L, cfg.K, cfg.gamma,
                    density, T, dT, m, alpha
                )
                if os.path.exists(sp + ".npz"):
                    results, _ = load_condition(sp)
                    data[T][dT][m] = results
                else:
                    missing_ms.append(m)

            found_ms = sorted(data[T][dT].keys())
            if found_ms:
                print(f"density={density:.4f}  T={T} dT={dT}: m={found_ms}")
            else:
                print(f"WARNING: no data found for density={density:.4f}, T={T}, dT={dT}")

            if missing_ms and found_ms:
                print(f"WARNING: missing m values for density={density:.4f}, T={T}, dT={dT}: {missing_ms}")

    return data, task_maps, alpha_maps


# ============================================================
# 5. FIGURE
# ============================================================

def fig2_m_sweep(
    data: Dict,
    task_maps: Dict,
    data_cfg: DataConfig,
    fig_cfg: FigureConfig,
    save_path: Optional[str] = None,
):
    """
    Flexible endpoint pheno_dist_norm vs m figure.

    One figure per density.
    One panel per T.
    One line per dT.
    """
    n_cols = len(data_cfg.T_values)
    fig_w = fig_cfg.figsize_per_col[0] * n_cols
    fig_h = fig_cfg.figsize_per_col[1]
    fig, axes = plt.subplots(1, n_cols, figsize=(fig_w, fig_h))

    if n_cols == 1:
        axes = [axes]

    fig.subplots_adjust(
        wspace=0.15,
        left=0.10, right=0.96, top=0.88, bottom=0.14
    )

    dt_color_map = make_dt_color_map(data_cfg.task_divs)

    target_dists = {
        T: {
            dT: mean_pairwise_distance(task_maps[T][dT])
            if (T in task_maps and dT in task_maps[T]) else np.nan
            for dT in data_cfg.task_divs
        }
        for T in data_cfg.T_values
    }

    for col, T in enumerate(data_cfg.T_values):
        ax = axes[col]
        ax.set_box_aspect(1)
        ax.set_title(f"T/K = {T / data_cfg.K:.1f}", fontsize=12, pad=6)

        all_ms_for_ticks = set()

        for dT in data_cfg.task_divs:
            if T not in data or dT not in data[T] or not data[T][dT]:
                continue

            td = target_dists[T][dT]
            if not (np.isfinite(td) and td > 0):
                print(f"WARNING: invalid normalization target for T={T}, dT={dT}")
                continue

            means = []
            sems = []
            ms_plot = []

            for m in sorted(data[T][dT].keys()):
                reps = data[T][dT][m]
                vals = endpoint_pheno_dist_norm(reps, td)
                if len(vals) == 0:
                    continue

                mean_val, sem_val = mean_and_sem(vals)
                means.append(mean_val)
                sems.append(sem_val)
                ms_plot.append(m)
                all_ms_for_ticks.add(m)

            if not ms_plot:
                continue

            ms_plot = np.array(ms_plot, dtype=float)
            means = np.array(means, dtype=float)
            sems = np.array(sems, dtype=float)

            color = dt_color_map[dT]
            label = make_dt_label(dT)

            ax.plot(
                ms_plot, means, "-o",
                color=color,
                label=label,
                lw=fig_cfg.line_width,
                ms=fig_cfg.marker_size,
                zorder=3
            )
            ax.fill_between(
                ms_plot,
                means - sems,
                means + sems,
                color=color,
                alpha=fig_cfg.band_alpha,
                zorder=2
            )

        ax.axhline(
            fig_cfg.reference_line_y,
            color="gray", ls="--", lw=0.8, alpha=0.5
        )

        ax.set_xlabel("m (tasks jointly evaluated)")
        ax.set_ylim(0, None)

        if all_ms_for_ticks:
            ax.set_xticks(sorted(all_ms_for_ticks))

        if col == 0:
            ax.set_ylabel("Realized differentiation\n(pheno dist / task separation)")
        else:
            ax.tick_params(labelleft=False)

    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles, labels,
            loc="upper center",
            ncol=min(len(labels), max(1, len(data_cfg.task_divs))),
            frameon=False,
            bbox_to_anchor=(0.5, 0.98),
            fontsize=10
        )

    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
        print(f"Saved: {save_path}")

    return fig


# ============================================================
# 6. NUMERICAL SUMMARY
# ============================================================

def print_summary(
    data: Dict,
    task_maps: Dict,
    data_cfg: DataConfig,
):
    """
    Print the exact values plotted:
    mean ± SEM of endpoint pheno_dist_norm for every (T/K, dT, m).
    """
    target_dists = {
        T: {
            dT: mean_pairwise_distance(task_maps[T][dT])
            if (T in task_maps and dT in task_maps[T]) else np.nan
            for dT in data_cfg.task_divs
        }
        for T in data_cfg.T_values
    }

    print("\n" + "=" * 68)
    print("NUMERICAL SUMMARY — endpoint pheno_dist_norm")
    print("=" * 68)
    print(f"{'T/K':>5}  {'dT':>5}  {'m':>3}  {'mean':>8}  {'sem':>8}  {'n_reps':>6}")
    print("-" * 68)

    for T in data_cfg.T_values:
        for dT in data_cfg.task_divs:
            if T not in data or dT not in data[T] or not data[T][dT]:
                continue

            td = target_dists[T][dT]
            if not (np.isfinite(td) and td > 0):
                print(f"WARNING: invalid normalization target for T={T}, dT={dT}")
                continue

            for m in sorted(data[T][dT].keys()):
                reps = data[T][dT][m]
                vals = endpoint_pheno_dist_norm(reps, td)
                if len(vals) == 0:
                    continue

                mean_val, sem_val = mean_and_sem(vals)
                print(
                    f"{T / data_cfg.K:>5.1f}  {dT:>5.1f}  {m:>3}  "
                    f"{mean_val:>8.4f}  {sem_val:>8.4f}  {len(vals):>6}"
                )


# ============================================================
# 7. CLI
# ============================================================

def parse_args():
    defaults = DataConfig()
    p = argparse.ArgumentParser()

    p.add_argument("--cache_dir", type=str, default=defaults.cache_dir)
    p.add_argument("--save_dir", type=str, default="figures_out")
    p.add_argument("--filename", type=str, default="fig2_m_sweep")
    p.add_argument("--fmt", type=str, default="pdf")

    p.add_argument("--L", type=int, default=defaults.L)
    p.add_argument("--K", type=int, default=defaults.K)
    p.add_argument("--gamma", type=float, default=defaults.gamma)

    p.add_argument("--densities", type=float, nargs="+",
                   default=defaults.densities, dest="densities")
    p.add_argument("--T", type=int, nargs="+",
                   default=defaults.T_values, dest="T_values")
    p.add_argument("--dT", type=float, nargs="+",
                   default=defaults.task_divs, dest="task_divs")

    p.add_argument("--no_show", action="store_true")
    p.add_argument("--no_summary", action="store_true")

    return p.parse_args()


# ============================================================
# 8. MAIN
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
        fig = fig2_m_sweep(
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