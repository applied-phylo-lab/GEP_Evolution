#!/usr/bin/env python3
"""
fig2_m_sweep.py
===============
Flexible plotting for cached simulation outputs.

Default behavior:
- One figure per density
- 3 rows x N columns, one column per T/K ratio
  Row 1: endpoint pheno_dist / mean pairwise task-target distance
  Row 2: endpoint usage cosine dissimilarity
  Row 3: endpoint Jaccard distance between program gene sets
- x-axis: m (simultaneity, 1 to T)
- lines: one per dT, grayscale color ramp (light = low dT, black = high dT)
- shaded band: mean ± SEM across replicates
- grayscale colorbar for dT (not legend)
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
    densities: List[float] = field(default_factory=lambda: [0.25])
    T_values: List[int] = field(default_factory=lambda: [2, 3, 4, 6, 8])
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
    figsize_per_col: Tuple[float, float] = (4.0, 3.8)
    line_width: float = 1.8
    marker_size: float = 5.0
    reference_line_y: float = 1.0


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
    """Endpoint normalized phenotype distance for each replicate."""
    if not (np.isfinite(target_dist) and target_dist > 0):
        return np.array([], dtype=float)
    vals = []
    for r in reps:
        n = r["n_actual_subs"]
        if n <= 0:
            continue
        vals.append(float(r["pheno_dist"][n - 1]) / target_dist)
    return np.array(vals, dtype=float)


def cosine_dissim_usage(usage: np.ndarray) -> float:
    """
    Mean pairwise cosine dissimilarity across task usage vectors.
    usage: (T, K), one row per task.
    """
    A = np.array(usage, dtype=float)
    norms = np.linalg.norm(A, axis=1, keepdims=True)
    norms = np.where(norms < 1e-12, 1.0, norms)
    A = A / norms
    vals = []
    for i, j in combinations(range(A.shape[0]), 2):
        vals.append(1.0 - float(np.dot(A[i], A[j])))
    return float(np.mean(vals)) if vals else 0.0


def jaccard_dist(genome: np.ndarray) -> float:
    """Mean pairwise Jaccard distance between program gene sets. genome: (L, K)."""
    K = genome.shape[1]
    if K < 2:
        return 0.0
    dists = []
    for j1, j2 in combinations(range(K), 2):
        g1 = genome[:, j1].astype(bool)
        g2 = genome[:, j2].astype(bool)
        u = np.sum(g1 | g2)
        dists.append(0.0 if u == 0 else 1.0 - np.sum(g1 & g2) / u)
    return float(np.mean(dists))


def endpoint_usage_cosine(reps: List[Dict]) -> np.ndarray:
    """Endpoint usage cosine dissimilarity for each replicate."""
    vals = []
    for r in reps:
        step = max(r["snapshots"].keys())
        usage = r["snapshots"][step]["usage"]
        vals.append(cosine_dissim_usage(usage))
    return np.array(vals, dtype=float)


def endpoint_jaccard(reps: List[Dict]) -> np.ndarray:
    """Endpoint Jaccard distance for each replicate."""
    vals = []
    for r in reps:
        step = max(r["snapshots"].keys())
        genome = r["snapshots"][step]["genome"]
        vals.append(jaccard_dist(genome))
    return np.array(vals, dtype=float)


def mean_and_std(vals: np.ndarray) -> Tuple[float, float]:
    """Return mean and std with safe behavior for small n."""
    if len(vals) == 0:
        return np.nan, np.nan
    if len(vals) == 1:
        return float(vals.mean()), 0.0
    return float(vals.mean()), float(vals.std(ddof=1))


# ============================================================
# 3. COLOR HELPERS
# ============================================================

def make_dt_gray_map(task_divs: List[float]) -> Dict[float, Tuple[float, float, float]]:
    """
    Map each dT to a grayscale RGB tuple.
    Linearly-spaced gray levels between 0.72 (lightest) and 0.0 (black).
    """
    vals = sorted(task_divs)
    if len(vals) == 1:
        return {vals[0]: (0.0, 0.0, 0.0)}
    gray_light = 0.72
    gray_dark = 0.0
    result = {}
    for i, dT in enumerate(vals):
        g = gray_light - (gray_light - gray_dark) * i / (len(vals) - 1)
        g = max(g, 0.0)
        result[dT] = (g, g, g)
    return result


def add_grayscale_colorbar(fig, task_divs: List[float]):
    """Single grayscale colorbar for dT, centered below the panels."""
    dt_min = min(task_divs)
    dt_max = max(task_divs)
    dt_norm = mpl.colors.Normalize(vmin=dt_min, vmax=dt_max)

    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "dt_gray", [(0.72, 0.72, 0.72), (0.0, 0.0, 0.0)]
    )

    cb_w = 0.50
    cb_x = 0.5 - cb_w / 2 + 0.02
    cb_y = 0.02
    cb_h = 0.018

    ax_cb = fig.add_axes([cb_x, cb_y, cb_w, cb_h])
    cb = mpl.colorbar.ColorbarBase(
        ax_cb, cmap=cmap, norm=dt_norm,
        orientation="horizontal", ticks=task_divs
    )
    cb.set_ticklabels([str(dT) for dT in task_divs])
    cb.set_label(r"$\Delta T$", fontsize=10)
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
            if missing_ms and found_ms:
                print(f"  missing m: {missing_ms}")

    return data, task_maps, alpha_maps


# ============================================================
# 5. FIGURE
# ============================================================

ROW_SPECS = [
    {
        "key": "pheno_dist_norm",
        "ylabel": r"Mean pairwise $\bar{\Delta Z} / \bar{\Delta T}$",
        "ylim": (0, None),
        "refline": 1.0,
    },
    {
        "key": "usage_cosine",
        "ylabel": r"Mean pairwise usage $\bar{\cos}$ dissimilarity",
        "ylim": (0, 1.05),
        "refline": None,
    },
    {
        "key": "jaccard",
        "ylabel": r"Mean pairwise Jaccard distance $\bar{J}$",
        "ylim": (0, 1.05),
        "refline": None,
    },
]


def _get_endpoint_vals(reps, row_key, target_dist=None):
    """Dispatch to the right endpoint extractor per row."""
    if row_key == "pheno_dist_norm":
        return endpoint_pheno_dist_norm(reps, target_dist)
    elif row_key == "usage_cosine":
        return endpoint_usage_cosine(reps)
    elif row_key == "jaccard":
        return endpoint_jaccard(reps)
    else:
        raise ValueError(f"Unknown row key: {row_key}")


def fig2_m_sweep(
    data: Dict,
    task_maps: Dict,
    data_cfg: DataConfig,
    fig_cfg: FigureConfig,
    save_path: Optional[str] = None,
):
    """
    3-row x N-col figure.
    Rows: pheno_dist_norm, usage cosine dissimilarity, Jaccard distance.
    Columns: one per T value.
    Lines: one per dT, grayscale.
    """
    n_rows = len(ROW_SPECS)
    n_cols = len(data_cfg.T_values)

    fig_w = fig_cfg.figsize_per_col[0] * n_cols
    fig_h = fig_cfg.figsize_per_col[1] * n_rows

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h))
    if n_cols == 1:
        axes = axes[:, np.newaxis]

    fig.subplots_adjust(
        wspace=0.18, hspace=0.32,
        left=0.10, right=0.96, top=0.93, bottom=0.08
    )

    dt_color_map = make_dt_gray_map(data_cfg.task_divs)

    # Precompute target distances
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

    for row, spec in enumerate(ROW_SPECS):
        for col, T in enumerate(data_cfg.T_values):
            ax = axes[row, col]
            ax.set_box_aspect(1)

            # Panel label
            if pi < len(panel_labels):
                ax.text(-0.15, 1.08, panel_labels[pi],
                        transform=ax.transAxes, fontsize=13, fontweight="bold",
                        va="top", ha="left")
            pi += 1

            # Column title (top row only)
            if row == 0:
                ax.set_title(f"T/K = {T / data_cfg.K:.1f}", fontsize=12, pad=6)

            all_ms_for_ticks = set()

            for dT in data_cfg.task_divs:
                if T not in data or dT not in data[T] or not data[T][dT]:
                    continue

                td = target_dists[T].get(dT, np.nan)

                means = []
                stds = []
                ms_plot = []

                for m in sorted(data[T][dT].keys()):
                    reps = data[T][dT][m]
                    vals = _get_endpoint_vals(reps, spec["key"], target_dist=td)

                    if len(vals) == 0:
                        continue

                    mean_val, std_val = mean_and_std(vals)
                    means.append(mean_val)
                    stds.append(std_val)
                    ms_plot.append(m)
                    all_ms_for_ticks.add(m)

                if not ms_plot:
                    continue

                ms_plot = np.array(ms_plot, dtype=float)
                means = np.array(means, dtype=float)
                stds = np.array(stds, dtype=float)

                color = dt_color_map[dT]

                ax.errorbar(
                    ms_plot, means,
                    yerr=stds,
                    fmt="-o",
                    color=color,
                    lw=fig_cfg.line_width,
                    ms=fig_cfg.marker_size,
                    capsize=2,
                    capthick=1.0,
                    elinewidth=1.0,
                    zorder=3,
                )

            # Reference line
            if spec["refline"] is not None:
                ax.axhline(spec["refline"], color="gray", ls="--", lw=0.8, alpha=0.5)

            # Axis formatting
            if spec["ylim"] is not None:
                ax.set_ylim(*spec["ylim"])

            if all_ms_for_ticks:
                ax.set_xticks(sorted(all_ms_for_ticks))

            # x-label only on bottom row
            if row == n_rows - 1:
                ax.set_xlabel("m (simultaneity)")
            else:
                ax.tick_params(labelbottom=False)

            # y-label only on leftmost column
            if col == 0:
                ax.set_ylabel(spec["ylabel"])
            else:
                ax.tick_params(labelleft=False)

    # Grayscale colorbar for dT
    add_grayscale_colorbar(fig, data_cfg.task_divs)

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
    Print exact values for all three metrics:
    pheno_dist_norm, usage cosine dissim, Jaccard distance.
    """
    target_dists = {
        T: {
            dT: mean_pairwise_distance(task_maps[T][dT])
            if (T in task_maps and dT in task_maps[T]) else np.nan
            for dT in data_cfg.task_divs
        }
        for T in data_cfg.T_values
    }

    print("\n" + "=" * 90)
    print("NUMERICAL SUMMARY — endpoint metrics")
    print("=" * 90)
    header = (
        f"{'T/K':>5}  {'dT':>5}  {'m':>3}  "
        f"{'pheno':>8}  {'p_sem':>7}  "
        f"{'cosine':>8}  {'c_sem':>7}  "
        f"{'jaccard':>8}  {'j_sem':>7}  "
        f"{'n':>4}"
    )
    print(header)
    print("-" * 90)

    for T in data_cfg.T_values:
        for dT in data_cfg.task_divs:
            if T not in data or dT not in data[T] or not data[T][dT]:
                continue

            td = target_dists[T].get(dT, np.nan)

            for m in sorted(data[T][dT].keys()):
                reps = data[T][dT][m]

                pv = endpoint_pheno_dist_norm(reps, td)
                cv = endpoint_usage_cosine(reps)
                jv = endpoint_jaccard(reps)

                pm, ps = mean_and_std(pv)
                cm, cs = mean_and_std(cv)
                jm, js = mean_and_std(jv)

                n = max(len(pv), len(cv), len(jv))

                print(
                    f"{T / data_cfg.K:>5.1f}  {dT:>5.1f}  {m:>3}  "
                    f"{pm:>8.4f}  {ps:>7.4f}  "
                    f"{cm:>8.4f}  {cs:>7.4f}  "
                    f"{jm:>8.4f}  {js:>7.4f}  "
                    f"{n:>4}"
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
            plt.close(fig)#!/usr/bin/env python3
"""
fig2_m_sweep.py
===============
Flexible plotting for cached simulation outputs.

Default behavior:
- One figure per density
- 3 rows x N columns, one column per T/K ratio
  Row 1: endpoint pheno_dist / mean pairwise task-target distance
  Row 2: endpoint usage cosine dissimilarity
  Row 3: endpoint Jaccard distance between program gene sets
- x-axis: m (simultaneity, 1 to T)
- lines: one per dT, grayscale color ramp (light = low dT, black = high dT)
- shaded band: mean ± SEM across replicates
- grayscale colorbar for dT (not legend)
"""

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


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
    filename: str = "fig2_m_sweep"
    fmt: str = "pdf"
    show: bool = True
    print_summary: bool = True


@dataclass
class FigureConfig:
    figsize_per_col: Tuple[float, float] = (4.0, 3.8)
    line_width: float = 1.8
    marker_size: float = 5.0
    reference_line_y: float = 1.0




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
            if missing_ms and found_ms:
                print(f"  missing m: {missing_ms}")

    return data, task_maps, alpha_maps


# ============================================================
# 5. FIGURE
# ============================================================

ROW_SPECS = [
    {
        "key": "pheno_dist_norm",
        "ylabel": r"Mean pairwise $\bar{\Delta Z} / \bar{\Delta T}$",
        "ylim": (0, None),
        "refline": 1.0,
    },
    {
        "key": "usage_cosine",
        "ylabel": r"Mean pairwise usage $\bar{\cos}$ dissimilarity",
        "ylim": (0, 1.05),
        "refline": None,
    },
    {
        "key": "jaccard",
        "ylabel": r"Mean pairwise Jaccard distance $\bar{J}$",
        "ylim": (0, 1.05),
        "refline": None,
    },
]




def fig2_m_sweep(
    data: Dict,
    task_maps: Dict,
    data_cfg: DataConfig,
    fig_cfg: FigureConfig,
    save_path: Optional[str] = None,
):
    """
    3-row x N-col figure.
    Rows: pheno_dist_norm, usage cosine dissimilarity, Jaccard distance.
    Columns: one per T value.
    Lines: one per dT, grayscale.
    """
    n_rows = len(ROW_SPECS)
    n_cols = len(data_cfg.T_values)

    fig_w = fig_cfg.figsize_per_col[0] * n_cols
    fig_h = fig_cfg.figsize_per_col[1] * n_rows

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h))
    if n_cols == 1:
        axes = axes[:, np.newaxis]

    fig.subplots_adjust(
        wspace=0.18, hspace=0.32,
        left=0.10, right=0.96, top=0.93, bottom=0.08
    )

    dt_color_map = make_dt_gray_map(data_cfg.task_divs)

    # Precompute target distances
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

    for row, spec in enumerate(ROW_SPECS):
        for col, T in enumerate(data_cfg.T_values):
            ax = axes[row, col]
            ax.set_box_aspect(1)

            # Panel label
            if pi < len(panel_labels):
                ax.text(-0.15, 1.08, panel_labels[pi],
                        transform=ax.transAxes, fontsize=13, fontweight="bold",
                        va="top", ha="left")
            pi += 1

            # Column title (top row only)
            if row == 0:
                ax.set_title(f"T/K = {T / data_cfg.K:.1f}", fontsize=12, pad=6)

            all_ms_for_ticks = set()

            for dT in data_cfg.task_divs:
                if T not in data or dT not in data[T] or not data[T][dT]:
                    continue

                td = target_dists[T].get(dT, np.nan)

                means = []
                stds = []
                ms_plot = []

                for m in sorted(data[T][dT].keys()):
                    reps = data[T][dT][m]
                    vals = _get_endpoint_vals(reps, spec["key"], target_dist=td)

                    if len(vals) == 0:
                        continue

                    mean_val, std_val = mean_and_std(vals)
                    means.append(mean_val)
                    stds.append(std_val)
                    ms_plot.append(m)
                    all_ms_for_ticks.add(m)

                if not ms_plot:
                    continue

                ms_plot = np.array(ms_plot, dtype=float)
                means = np.array(means, dtype=float)
                stds = np.array(stds, dtype=float)

                color = dt_color_map[dT]

                ax.errorbar(
                    ms_plot, means,
                    yerr=stds,
                    fmt="-o",
                    color=color,
                    lw=fig_cfg.line_width,
                    ms=fig_cfg.marker_size,
                    capsize=2,
                    capthick=1.0,
                    elinewidth=1.0,
                    zorder=3,
                )

            # Reference line
            if spec["refline"] is not None:
                ax.axhline(spec["refline"], color="gray", ls="--", lw=0.8, alpha=0.5)

            # Axis formatting
            if spec["ylim"] is not None:
                ax.set_ylim(*spec["ylim"])

            if all_ms_for_ticks:
                ax.set_xticks(sorted(all_ms_for_ticks))

            # x-label only on bottom row
            if row == n_rows - 1:
                ax.set_xlabel("m (simultaneity)")
            else:
                ax.tick_params(labelbottom=False)

            # y-label only on leftmost column
            if col == 0:
                ax.set_ylabel(spec["ylabel"])
            else:
                ax.tick_params(labelleft=False)

    # Grayscale colorbar for dT
    add_grayscale_colorbar(fig, data_cfg.task_divs)

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
    Print exact values for all three metrics:
    pheno_dist_norm, usage cosine dissim, Jaccard distance.
    """
    target_dists = {
        T: {
            dT: mean_pairwise_distance(task_maps[T][dT])
            if (T in task_maps and dT in task_maps[T]) else np.nan
            for dT in data_cfg.task_divs
        }
        for T in data_cfg.T_values
    }

    print("\n" + "=" * 90)
    print("NUMERICAL SUMMARY — endpoint metrics")
    print("=" * 90)
    header = (
        f"{'T/K':>5}  {'dT':>5}  {'m':>3}  "
        f"{'pheno':>8}  {'p_sem':>7}  "
        f"{'cosine':>8}  {'c_sem':>7}  "
        f"{'jaccard':>8}  {'j_sem':>7}  "
        f"{'n':>4}"
    )
    print(header)
    print("-" * 90)

    for T in data_cfg.T_values:
        for dT in data_cfg.task_divs:
            if T not in data or dT not in data[T] or not data[T][dT]:
                continue

            td = target_dists[T].get(dT, np.nan)

            for m in sorted(data[T][dT].keys()):
                reps = data[T][dT][m]

                pv = endpoint_pheno_dist_norm(reps, td)
                cv = endpoint_usage_cosine(reps)
                jv = endpoint_jaccard(reps)

                pm, ps = mean_and_std(pv)
                cm, cs = mean_and_std(cv)
                jm, js = mean_and_std(jv)

                n = max(len(pv), len(cv), len(jv))

                print(
                    f"{T / data_cfg.K:>5.1f}  {dT:>5.1f}  {m:>3}  "
                    f"{pm:>8.4f}  {ps:>7.4f}  "
                    f"{cm:>8.4f}  {cs:>7.4f}  "
                    f"{jm:>8.4f}  {js:>7.4f}  "
                    f"{n:>4}"
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