#!/usr/bin/env python3
"""
fig_2_modularity.py
===================

Companion to fig_2.py. Same column structure (m=1 vs m=T) and same
visual style, but plots two complementary notions of modularity at the
fixed substitution cutoff.

Row 1 — Deployment / regulatory modularity:
    Mean pairwise cosine dissimilarity (1 - cos) across the program-usage
    vectors that the genome deploys for each task. Reads as "how
    differently does the genome use its programs across tasks?"
    0 -> identical deployment across tasks (no regulatory DOL).
    1 -> orthogonal deployment (full regulatory DOL).

Row 2 — Mutational / variational modularity (Tikhonov-style M):
    For each task pair (a, b), compute the cloud of single-bit mutation
    effects (dF_a, dF_b) over all L*K mutations. In polar (rho, phi),
        M_ab = 1 - sum(rho * |sin(2 phi)|) / sum(rho).
    Average over task pairs. Reads as "are mutational effects aligned
    with task axes, or mixed across tasks?"
    1 -> mutations are task-orthogonal (modular variation).
    0 -> mutations affect tasks equally (non-modular variation).

Layout: 2 rows x 2 cols.
    cols : m = 1 (sequential), m = T (simultaneous)
    rows : usage cosine dissimilarity, mutational modularity
    x-axis: T / K
    lines : one per dT, viridis_r color ramp (matches fig_2.py)

Mutational-modularity computation reuses helpers from
figS_tradeoff_modularity to ensure the M definition is identical.
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

# Reuse the exact M definition used in figS_tradeoff_modularity to keep the
# two figures consistent.
from figS_tradeoff_modularity import mean_pairwise_tradeoff_modularity


# ============================================================
# 0. MPL STYLE  (matches fig_2.py)
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
    cutoff: int = 200


@dataclass
class OutputConfig:
    save_dir: str = "figures_out"
    filename: str = "fig_2_modularity"
    fmt: str = "pdf"
    show: bool = True
    print_summary: bool = True


@dataclass
class GroupSpec:
    label: str
    m_selector: MSelector


@dataclass
class FigureConfig:
    figsize: Tuple[float, float] = (9.0, 8.0)
    groups: List[GroupSpec] = field(default_factory=lambda: [
        GroupSpec(label="m = 1", m_selector="min"),
        GroupSpec(label="m = T", m_selector="T"),
    ])


# ============================================================
# 2. PER-METRIC HELPERS
# ============================================================

def usage_cosine_dissim(usage: np.ndarray) -> float:
    """
    Mean pairwise cosine dissimilarity across task usage vectors.
    usage : (T, K) -- one row per task.
    0 -> identical deployment across tasks; 1 -> orthogonal deployment.
    """
    A = np.array(usage, dtype=float)
    norms = np.linalg.norm(A, axis=1, keepdims=True)
    norms = np.where(norms < 1e-12, 1.0, norms)
    A = A / norms
    pairs = list(combinations(range(A.shape[0]), 2))
    if not pairs:
        return np.nan
    return float(np.mean([1.0 - np.dot(A[i], A[j]) for i, j in pairs]))


def _snapshot_at_cutoff(r: Dict, cutoff: int) -> Dict:
    """Snapshot at step closest to but not exceeding cutoff-1."""
    snap_steps = sorted(r["snapshots"].keys())
    eligible = [s for s in snap_steps if s <= cutoff - 1]
    step = eligible[-1] if eligible else snap_steps[0]
    return r["snapshots"][step]


def usage_cosine_at_step(reps: List[Dict], cutoff: int) -> np.ndarray:
    """Per-replicate usage cosine dissimilarity at snapshot <= cutoff."""
    return np.array(
        [usage_cosine_dissim(_snapshot_at_cutoff(r, cutoff)["usage"])
         for r in reps],
        dtype=float,
    )


def mutational_M_at_step(
    reps: List[Dict],
    tasks: np.ndarray,
    gamma: float,
    cutoff: int,
) -> np.ndarray:
    """
    Per-replicate Tikhonov-style mean-pairwise mutational modularity M
    at snapshot <= cutoff. Reuses mean_pairwise_tradeoff_modularity from
    figS_tradeoff_modularity.

    Returns array of M values (one per replicate).
    """
    vals = []
    for r in reps:
        snap = _snapshot_at_cutoff(r, cutoff)
        # mean_pairwise_tradeoff_modularity returns (M, chi); we only want M.
        M, _ = mean_pairwise_tradeoff_modularity(snap["genome"], tasks, gamma)
        vals.append(M)
    return np.array(vals, dtype=float)


def mean_and_std(vals: np.ndarray) -> Tuple[float, float]:
    finite = vals[np.isfinite(vals)] if len(vals) > 0 else vals
    if len(finite) == 0:
        return np.nan, np.nan
    if len(finite) == 1:
        return float(finite.mean()), 0.0
    return float(finite.mean()), float(finite.std(ddof=1))


# ============================================================
# 3. COLOR HELPERS  (matches fig_2.py — viridis_r)
# ============================================================

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
    ax_cb = fig.add_axes([0.5 - cb_w / 2 + 0.02, 0.02, cb_w, 0.018])
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


# ============================================================
# 4. DATA LOADING  (mirrors fig_2.py)
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
# 5. SHARED PLOT LOOP  (mirrors fig_2.py _plot_lines)
# ============================================================

def _plot_lines(ax, data, data_cfg, fig_cfg, group, dt_color_map, get_vals_fn):
    """Plot errorbar lines over T/K for all dT into ax.

    get_vals_fn(reps, T, dT) -> 1D array of per-replicate metric values.
    """
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


# ============================================================
# 6. FIGURE
# ============================================================

ROW_SPECS = [
    {
        "key": "usage_cosine",
        "ylabel": r"Deployment modularity, $1 - \overline{\cos}$",
        "ylim": (0, 1.05),
    },
    {
        "key": "mutational_M",
        "ylabel": r"Mutational modularity, $\bar{M}$",
        "ylim": (0, 1.05),
    },
]


def _vals_for_row(reps, row_key, cutoff, tasks=None, gamma=None):
    if row_key == "usage_cosine":
        return usage_cosine_at_step(reps, cutoff)
    if row_key == "mutational_M":
        if tasks is None or gamma is None:
            raise ValueError("mutational_M requires tasks and gamma")
        return mutational_M_at_step(reps, tasks, gamma, cutoff)
    raise ValueError(f"Unknown row key: {row_key}")


def fig_modularity(data, task_maps, data_cfg, fig_cfg, cutoff, save_path=None):
    """2 x 2 figure: rows = usage cosine, mutational M;  cols = m=1, m=T."""
    n_groups = len(fig_cfg.groups)
    n_rows = len(ROW_SPECS)
    fig, axes = plt.subplots(n_rows, n_groups, figsize=fig_cfg.figsize)
    if n_groups == 1:
        axes = axes[:, np.newaxis]
    fig.subplots_adjust(wspace=0.25, hspace=0.30,
                        left=0.12, right=0.94, top=0.92, bottom=0.12)

    dt_color_map = make_dt_viridis_map(data_cfg.task_divs)
    panel_labels = "ABCD"

    for row, spec in enumerate(ROW_SPECS):
        for col, group in enumerate(fig_cfg.groups):
            ax = axes[row, col]
            ax.set_box_aspect(1)
            ax.text(-0.18, 1.08, panel_labels[row * n_groups + col],
                    transform=ax.transAxes, fontsize=14, fontweight="bold",
                    va="top", ha="left")
            if row == 0:
                ax.set_title(f"{group.label}  (sub $\\leq$ {cutoff})",
                             fontsize=12, pad=6)

            def make_getter(row_key=spec["key"]):
                def fn(reps, T, dT):
                    tasks = task_maps.get(T, {}).get(dT, None)
                    return _vals_for_row(
                        reps, row_key, cutoff,
                        tasks=tasks, gamma=data_cfg.gamma,
                    )
                return fn

            _plot_lines(
                ax, data, data_cfg, fig_cfg, group, dt_color_map,
                make_getter(),
            )

            if spec["ylim"] is not None:
                ax.set_ylim(*spec["ylim"])
            if col == 0:
                ax.set_ylabel(spec["ylabel"])
            else:
                ax.tick_params(labelleft=False)
            if row == n_rows - 1:
                ax.set_xlabel("T / K")
            else:
                ax.tick_params(labelbottom=False)

    _add_viridis_colorbar(fig, data_cfg.task_divs)
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
        print(f"Saved: {save_path}")
    return fig


# ============================================================
# 7. NUMERICAL SUMMARY
# ============================================================

def print_summary(data, task_maps, data_cfg, fig_cfg, cutoff):
    print(f"\n{'=' * 92}")
    print(f"SUMMARY — modularity at sub ≤ {cutoff}")
    print("=" * 92)
    header = (f"{'T/K':>6}  {'dT':>5}  {'group':>8}  {'m':>4}  "
              f"{'cos_d':>8}  {'cd_std':>7}  "
              f"{'mut_M':>8}  {'mM_std':>7}  {'n':>4}")
    print(header)
    print("-" * len(header))

    for T in data_cfg.T_values:
        for dT in data_cfg.task_divs:
            if dT not in data[T] or not data[T][dT]:
                continue
            available_m = sorted(data[T][dT].keys())
            tasks = task_maps.get(T, {}).get(dT, None)

            for group in fig_cfg.groups:
                m = resolve_m_selector(group.m_selector, available_m, T)
                if m is None or m not in data[T][dT]:
                    continue
                reps = data[T][dT][m]

                cv = usage_cosine_at_step(reps, cutoff)
                if tasks is not None:
                    Mv = mutational_M_at_step(reps, tasks, data_cfg.gamma, cutoff)
                else:
                    Mv = np.array([], dtype=float)
                cm, cs = mean_and_std(cv)
                Mm, Ms = mean_and_std(Mv)

                print(f"{T / data_cfg.K:>6.1f}  {dT:>5.1f}  {group.label:>8}  {m:>4}  "
                      f"{cm:>8.4f}  {cs:>7.4f}  "
                      f"{Mm:>8.4f}  {Ms:>7.4f}  {len(cv):>4}")


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

        save_path = os.path.join(
            out_cfg.save_dir,
            f"{out_cfg.filename}_sub{data_cfg.cutoff}"
            f"_gamma{data_cfg.gamma}_fr{data_cfg.fitness_r}"
            f"_density{density:.4f}.{out_cfg.fmt}"
        )

        fig = fig_modularity(data, task_maps, data_cfg, fig_cfg,
                             cutoff=data_cfg.cutoff, save_path=save_path)

        if out_cfg.print_summary:
            print_summary(data, task_maps, data_cfg, fig_cfg, cutoff=data_cfg.cutoff)

        if out_cfg.show:
            plt.show()
        else:
            plt.close(fig)