#!/usr/bin/env python3
"""
fig_diagnostics.py
==================

Diagnostic trajectories for identifying the early-time regime.

Four figures (one per metric): W, mean P, n_ben, s_max.
Layout: N_T rows x 2 cols (m=1 | m=T).
Trajectories truncated at `max_sub` substitutions.
Pass --max_sub 200 to inspect up to step 200.
"""

import argparse
import json
import os
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

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
    max_sub: int = 200


@dataclass
class OutputConfig:
    save_dir: str = "figures_out"
    filename: str = "fig_diagnostics"
    fmt: str = "pdf"
    show: bool = True


@dataclass
class FigureConfig:
    figsize_per_panel: Tuple[float, float] = (3.5, 3.0)
    line_width: float = 0.75
    alpha: float = 0.8
    downsample: int = 10


# ============================================================
# 2. HELPERS
# ============================================================

def make_dt_viridis_map(task_divs: List[float]) -> Dict[float, Tuple]:
    vals = sorted(task_divs)
    cmap = mpl.colormaps["viridis_r"]
    if len(vals) == 1:
        return {vals[0]: cmap(0.5)}
    return {dT: cmap(i / (len(vals) - 1)) for i, dT in enumerate(vals)}


def _add_colorbar_horizontal(fig, task_divs: List[float]):
    dt_norm = mpl.colors.Normalize(vmin=min(task_divs), vmax=max(task_divs))
    cmap = mpl.colormaps["viridis_r"]
    ax_cb = fig.add_axes([0.5 - 0.25 + 0.02, 0.02, 0.50, 0.018])
    cb = mpl.colorbar.ColorbarBase(
        ax_cb, cmap=cmap, norm=dt_norm,
        orientation="horizontal", ticks=task_divs
    )
    cb.set_ticklabels([str(dT) for dT in task_divs])
    cb.set_label(r"$\Delta T$", fontsize=10)
    cb.ax.tick_params(labelsize=8)


# ============================================================
# 3. DATA LOADING
# ============================================================

def load_all_for_density(cfg: DataConfig, density: float):
    data, task_maps, alpha_maps = {}, {}, {}

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
                sp = sim_cache_path(cfg.cache_dir, cfg.L, cfg.K, cfg.gamma,
                                    density, T, dT, m, alpha)
                if os.path.exists(sp + ".npz"):
                    results, _ = load_condition(sp)
                    data[T][dT][m] = results

        for dT in cfg.task_divs:
            if data[T].get(dT):
                print(f"density={density:.4f}  T={T} dT={dT}: m={sorted(data[T][dT].keys())}")

    return data, task_maps, alpha_maps


# ============================================================
# 4. TRAJECTORY EXTRACTION
# ============================================================

def _align_and_stats(
    trajs: List[np.ndarray], downsample: int, max_sub: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Truncate each trajectory to max_sub, NaN-pad to common length,
    downsample, return (x, mean, std).
    """
    if not trajs:
        return np.array([]), np.array([]), np.array([])

    trajs = [t[:max_sub] for t in trajs]
    max_n = max(len(t) for t in trajs)
    mat = np.full((len(trajs), max_n), np.nan)
    for i, t in enumerate(trajs):
        mat[i, :len(t)] = t

    x    = np.arange(max_n)[::downsample]
    mat  = mat[:, ::downsample]
    mean = np.nanmean(mat, axis=0)
    std  = np.nanstd(mat, axis=0, ddof=1)
    return x, mean, std


def _extract(reps: List[Dict], key: str, downsample: int, max_sub: int,
             transform: Optional[Callable] = None):
    trajs = []
    for r in reps:
        n = min(r["n_actual_subs"], max_sub)
        arr = np.array(r[key][:n], dtype=float)
        trajs.append(transform(arr) if transform else arr)
    return _align_and_stats(trajs, downsample, max_sub)


# ============================================================
# 5. METRIC SPECS
# ============================================================

METRIC_SPECS = [
    {
        "key":       "P",
        "transform": lambda P: P.mean(axis=1),   # (n, T) -> (n,)
        "ylabel":    r"Mean performance $\bar{P}$",
        "tag":       "P",
        "refline":   None,
        "ylim":      (0, 1.05),
    },
    {
        "key":       "n_ben",
        "transform": None,   # divided by L*K inside _make_diagnostic_fig
        "ylabel":    r"Fraction beneficial mutations $n_\mathrm{ben} / LK$",
        "tag":       "nben",
        "refline":   None,
        "ylim":      (0, None),
    },
]


# ============================================================
# 6. FIGURE BUILDER
# ============================================================

def _make_diagnostic_fig(
    data: Dict,
    data_cfg: DataConfig,
    fig_cfg: FigureConfig,
    spec: Dict,
    save_path: Optional[str] = None,
) -> plt.Figure:
    n_rows = len(data_cfg.T_values)
    fig_w = fig_cfg.figsize_per_panel[0] * 2
    fig_h = fig_cfg.figsize_per_panel[1] * n_rows

    fig, axes = plt.subplots(n_rows, 2, figsize=(fig_w, fig_h))
    if n_rows == 1:
        axes = axes[np.newaxis, :]

    fig.subplots_adjust(wspace=0.10, hspace=0.30,
                        left=0.14, right=0.96, top=0.95, bottom=0.08)

    dt_color_map = make_dt_viridis_map(data_cfg.task_divs)

    for row, T in enumerate(data_cfg.T_values):
        for col, (m, col_label) in enumerate([(1, "m = 1"), (T, "m = T")]):
            ax = axes[row, col]
            ax.set_box_aspect(1)

            if row == 0:
                ax.set_title(col_label, fontsize=12, pad=6)
            if col == 0:
                ax.set_ylabel(f"T/K = {T / data_cfg.K:.1f}\n" + spec["ylabel"],
                              fontsize=10)
            else:
                ax.tick_params(labelleft=False)
            if row == n_rows - 1:
                ax.set_xlabel("Substitution step")
            else:
                ax.tick_params(labelbottom=False)

            for dT in data_cfg.task_divs:
                if dT not in data[T] or m not in data[T][dT]:
                    continue
                color = dt_color_map[dT]
                transform = spec["transform"]
                if spec["key"] == "n_ben":
                    LK = data_cfg.L * data_cfg.K
                    transform = lambda arr: arr / LK
                x, mean, std = _extract(
                    data[T][dT][m], spec["key"],
                    fig_cfg.downsample, data_cfg.max_sub,
                    transform,
                )
                if len(x) == 0:
                    continue
                good = np.isfinite(mean)
                ax.errorbar(
                    x[good], mean[good], yerr=std[good],
                    fmt="-o", color=color,
                    lw=fig_cfg.line_width, ms=3,
                    markerfacecolor="none", markeredgecolor=color,
                    capsize=2, capthick=0.8, elinewidth=0.8,
                    alpha=fig_cfg.alpha,
                )

            if spec["refline"] is not None:
                ax.axhline(spec["refline"], color="gray", ls="--", lw=0.8, alpha=0.5)
            if spec["ylim"] is not None:
                ax.set_ylim(*spec["ylim"])

    _add_colorbar_horizontal(fig, data_cfg.task_divs)
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
        print(f"Saved: {save_path}")
    return fig



# ============================================================
# 6b. FIGURE BUILDER — rows=T/K, cols=m
# ============================================================

def _make_diagnostic_fig_m_sweep(
    data: Dict,
    data_cfg: DataConfig,
    fig_cfg: FigureConfig,
    spec: Dict,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Same metrics but rows=T/K, cols=all available m."""
    all_m = sorted({
        m for T in data_cfg.T_values
        for dT in data_cfg.task_divs
        if dT in data[T]
        for m in data[T][dT].keys()
    })
    if not all_m:
        print("No data found.")
        return None

    n_rows = len(data_cfg.T_values)
    n_cols = len(all_m)
    fig_w = fig_cfg.figsize_per_panel[0] * n_cols
    fig_h = fig_cfg.figsize_per_panel[1] * n_rows

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h))
    if n_rows == 1:
        axes = axes[np.newaxis, :]
    if n_cols == 1:
        axes = axes[:, np.newaxis]

    fig.subplots_adjust(wspace=0.10, hspace=0.30,
                        left=0.14, right=0.96, top=0.95, bottom=0.08)

    dt_color_map = make_dt_viridis_map(data_cfg.task_divs)

    for row, T in enumerate(data_cfg.T_values):
        for col, m in enumerate(all_m):
            ax = axes[row, col]

            # hide panels where m > T (no data possible)
            if m > T:
                ax.set_visible(False)
                continue

            ax.set_box_aspect(1)

            if row == 0:
                ax.set_title(f"m = {m}", fontsize=12, pad=6)
            if col == 0:
                ax.set_ylabel(f"T/K = {T / data_cfg.K:.1f}\n" + spec["ylabel"],
                              fontsize=10)
            else:
                ax.tick_params(labelleft=False)
            if row == n_rows - 1:
                ax.set_xlabel("Substitution step")
            else:
                ax.tick_params(labelbottom=False)

            for dT in data_cfg.task_divs:
                if dT not in data[T] or m not in data[T][dT]:
                    continue
                color = dt_color_map[dT]
                transform = spec["transform"]
                if spec["key"] == "n_ben":
                    LK = data_cfg.L * data_cfg.K
                    transform = lambda arr: arr / LK
                x, mean, std = _extract(
                    data[T][dT][m], spec["key"],
                    fig_cfg.downsample, data_cfg.max_sub,
                    transform,
                )
                if len(x) == 0:
                    continue
                good = np.isfinite(mean)
                ax.errorbar(
                    x[good], mean[good], yerr=std[good],
                    fmt="-o", color=color,
                    lw=fig_cfg.line_width, ms=3,
                    markerfacecolor="none", markeredgecolor=color,
                    capsize=2, capthick=0.8, elinewidth=0.8,
                    alpha=fig_cfg.alpha,
                )

            if spec["refline"] is not None:
                ax.axhline(spec["refline"], color="gray", ls="--", lw=0.8, alpha=0.5)
            if spec["ylim"] is not None:
                ax.set_ylim(*spec["ylim"])

    _add_colorbar_horizontal(fig, data_cfg.task_divs)
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
        print(f"Saved: {save_path}")
    return fig

# ============================================================
# 7. CLI + MAIN
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
    p.add_argument("--densities",  type=float, nargs="+", default=d.densities)
    p.add_argument("--T",          type=int,   nargs="+", default=d.T_values,  dest="T_values")
    p.add_argument("--dT",         type=float, nargs="+", default=d.task_divs, dest="task_divs")
    p.add_argument("--max_sub",    type=int,   default=d.max_sub)
    p.add_argument("--no_show",    action="store_true")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()

    data_cfg = DataConfig(
        cache_dir=args.cache_dir, L=args.L, K=args.K, gamma=args.gamma,
        densities=args.densities, T_values=args.T_values,
        task_divs=args.task_divs, max_sub=args.max_sub,
    )
    out_cfg = OutputConfig(
        save_dir=args.save_dir, filename=args.filename,
        fmt=args.fmt, show=not args.no_show,
    )
    fig_cfg = FigureConfig()
    os.makedirs(out_cfg.save_dir, exist_ok=True)

    for density in data_cfg.densities:
        print(f"\nLoading data for density={density:.4f}...")
        data, _, _ = load_all_for_density(data_cfg, density)

        for spec in METRIC_SPECS:
            # Figure A: rows=T/K, cols=m=1|m=T
            save_path = os.path.join(
                out_cfg.save_dir,
                f"{out_cfg.filename}_{spec['tag']}_sub{data_cfg.max_sub}_density{density:.4f}.{out_cfg.fmt}"
            )
            print(f"Plotting {spec['key']} (rows=T/K, cols=m=1|m=T)...")
            fig = _make_diagnostic_fig(
                data=data, data_cfg=data_cfg, fig_cfg=fig_cfg,
                spec=spec, save_path=save_path,
            )
            if out_cfg.show:
                plt.show()
            else:
                plt.close(fig)

            # Figure B: rows=T/K, cols=all m
            save_path_m = os.path.join(
                out_cfg.save_dir,
                f"{out_cfg.filename}_{spec['tag']}_msweep_sub{data_cfg.max_sub}_density{density:.4f}.{out_cfg.fmt}"
            )
            print(f"Plotting {spec['key']} (rows=T/K, cols=all m)...")
            fig_m = _make_diagnostic_fig_m_sweep(
                data=data, data_cfg=data_cfg, fig_cfg=fig_cfg,
                spec=spec, save_path=save_path_m,
            )
            if out_cfg.show:
                plt.show()
            elif fig_m is not None:
                plt.close(fig_m)