#!/usr/bin/env python3
"""
figS_tradeoff_modularity.py
===========================

Post hoc plotting of a Tikhonov-style tradeoff-modularity plane from
cached snapshot genomes only. No simulation reruns.

Now includes a thin cache layer for the computed plane coordinates.

Figure layout
-------------
- rows: m = 1 and m = T
- columns: one per T/K
- within each subplot:
    - faint viridis trajectories for all replicates of each dT condition
    - one colored mean trajectory per dT
    - colors follow the same viridis_r dT colormap used elsewhere

Per snapshot metric
-------------------
For each snapshot genome:
1. compute fitness/performance on all T tasks
2. enumerate all single-bit mutations
3. get dF matrix of shape (L*K, T)
4. for each task pair (a, b), compute:
      - Tikhonov-style mutational tradeoff chi_ab
      - Tikhonov-style mutational modularity M_ab
5. average over all task pairs:
      mean_pairwise_tradeoff
      mean_pairwise_modularity

These averaged coordinates define one point in the plane for that snapshot.

Plane-cache layer
-----------------
For each condition, save a compact cache:
- trajs: (n_reps, n_points, 2)
- snapshot_steps_per_rep: padded integer array
- params metadata in json

This avoids recomputing mutation clouds every time you redraw figures.
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

from simulate import (
    load_condition,
    load_task_map,
    sim_cache_path,
    task_cache_path,
    compute_performance,
)


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
    plane_cache_dir: str = "plane_cache"
    L: int = 100
    K: int = 4
    gamma: float = 1.0
    fitness_r: float = 0.0
    density: float = 0.25
    T_values: List[int] = field(default_factory=lambda: [2, 3, 4, 6, 8])
    task_divs: List[float] = field(default_factory=lambda: [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])


@dataclass
class FigureConfig:
    figsize_per_col: Tuple[float, float] = (4.2, 4.2)
    faint_alpha: float = 0.28
    faint_ms: float = 2.0
    mean_lw: float = 1.8
    arrow_lw: float = 1.8
    panel_label_fs: float = 13
    title_fs: float = 12
    carry_forward_to_n_points: int = 12


# ============================================================
# 2. COLOR HELPERS
# ============================================================

def make_dt_viridis_map(task_divs: List[float]) -> Dict[float, Tuple]:
    """Map each dT to viridis_r color, linearly spaced."""
    vals = sorted(task_divs)
    cmap = mpl.colormaps["viridis_r"]
    if len(vals) == 1:
        return {vals[0]: cmap(0.5)}
    return {dT: cmap(i / (len(vals) - 1)) for i, dT in enumerate(vals)}


def add_viridis_colorbar(fig, task_divs: List[float]):
    dt_norm = mpl.colors.Normalize(vmin=min(task_divs), vmax=max(task_divs))
    cmap = mpl.colormaps["viridis_r"]

    cb_w = 0.45
    cb_x = 0.5 - cb_w / 2
    cb_y = 0.03
    cb_h = 0.018

    ax_cb = fig.add_axes([cb_x, cb_y, cb_w, cb_h])
    cb = mpl.colorbar.ColorbarBase(
        ax_cb, cmap=cmap, norm=dt_norm,
        orientation="horizontal", ticks=task_divs
    )
    cb.set_ticklabels([str(x) for x in task_divs])
    cb.set_label(r"$\Delta T$", fontsize=10)
    cb.ax.tick_params(labelsize=8)


# ============================================================
# 3. CACHE LOADING
# ============================================================

def load_alpha_map(cache_dir: str, L: int, K: int, gamma: float,
                   fitness_r: float, T: int) -> Dict[float, float]:
    tpath = task_cache_path(cache_dir, L, K, gamma, fitness_r, T)
    with open(tpath + "_meta.json") as f:
        meta = json.load(f)
    return {float(k): v for k, v in meta["alpha_map"].items()}


def load_condition_if_exists(
    cache_dir: str,
    L: int,
    K: int,
    gamma: float,
    fitness_r: float,
    density: float,
    T: int,
    dT: float,
    m: int,
    alpha: float,
):
    sp = sim_cache_path(cache_dir, L, K, gamma, fitness_r, density, T, dT, m, alpha)
    if not os.path.exists(sp + ".npz"):
        return None, None
    return load_condition(sp)


# ============================================================
# 4. PLANE CACHE HELPERS
# ============================================================

def _plane_param_root(plane_cache_dir: str, L: int, K: int, gamma: float,
                      fitness_r: float) -> str:
    folder = os.path.join(plane_cache_dir, f"L{L}_K{K}_gamma{gamma}_fr{fitness_r}")
    os.makedirs(folder, exist_ok=True)
    return folder


def _plane_density_root(plane_cache_dir: str, L: int, K: int, gamma: float,
                        fitness_r: float, density: float) -> str:
    folder = os.path.join(
        _plane_param_root(plane_cache_dir, L, K, gamma, fitness_r),
        f"density{density:.4f}"
    )
    os.makedirs(folder, exist_ok=True)
    return folder


def plane_cache_path(
    plane_cache_dir: str,
    L: int,
    K: int,
    gamma: float,
    fitness_r: float,
    density: float,
    T: int,
    dT: float,
    m: int,
    alpha: float,
) -> str:
    return os.path.join(
        _plane_density_root(plane_cache_dir, L, K, gamma, fitness_r, density),
        f"plane_T{T}_dT{dT}_m{m}_alpha{alpha:.4f}"
    )


def save_plane_condition(base_path: str, trajs: np.ndarray, snapshot_steps_per_rep: List[List[int]], meta: Dict):
    """
    Save compact plane cache for one condition.

    trajs shape: (n_reps, n_points, 2)
    snapshot_steps_per_rep: ragged list, padded to rectangular integer array
    """
    n_reps = len(snapshot_steps_per_rep)
    max_len = max((len(x) for x in snapshot_steps_per_rep), default=0)
    padded_steps = -np.ones((n_reps, max_len), dtype=int)

    for i, steps in enumerate(snapshot_steps_per_rep):
        if len(steps) > 0:
            padded_steps[i, :len(steps)] = np.array(steps, dtype=int)

    np.savez_compressed(
        base_path + ".npz",
        trajs=trajs,
        snapshot_steps_per_rep=padded_steps,
    )

    with open(base_path + "_meta.json", "w") as f:
        json.dump(meta, f, indent=2)


def load_plane_condition(base_path: str) -> Tuple[np.ndarray, List[List[int]], Dict]:
    arr = np.load(base_path + ".npz", allow_pickle=False)
    trajs = arr["trajs"]
    padded_steps = arr["snapshot_steps_per_rep"]

    snapshot_steps_per_rep = []
    for row in padded_steps:
        snapshot_steps_per_rep.append([int(x) for x in row if int(x) >= 0])

    with open(base_path + "_meta.json") as f:
        meta = json.load(f)

    return trajs, snapshot_steps_per_rep, meta


# ============================================================
# 5. PAIRWISE TIKHONOV-STYLE METRICS
# ============================================================

def pairwise_tradeoff_for_two_tasks(dFa: np.ndarray, dFb: np.ndarray) -> float:
    ben = (dFa > 0) | (dFb > 0)
    if not np.any(ben):
        return np.nan

    xa = dFa[ben]
    xb = dFb[ben]

    denom = np.sqrt(np.sum(xa ** 2) * np.sum(xb ** 2))
    if denom <= 1e-15:
        return np.nan

    chi = -np.sum(xa * xb) / denom
    return float(np.clip(chi, -1.0, 1.0))


def pairwise_modularity_for_two_tasks(dFa: np.ndarray, dFb: np.ndarray) -> float:
    rho = np.sqrt(dFa ** 2 + dFb ** 2)
    valid = rho > 1e-15
    if not np.any(valid):
        return np.nan

    phi = np.arctan2(dFb[valid], dFa[valid])
    numer = np.sum(rho[valid] * np.abs(np.sin(2.0 * phi)))
    denom = np.sum(rho[valid])

    if denom <= 1e-15:
        return np.nan

    M = 1.0 - numer / denom
    return float(np.clip(M, 0.0, 1.0))


def mutation_effect_matrix(genome: np.ndarray, tasks: np.ndarray, gamma: float) -> np.ndarray:
    wt_info = compute_performance(genome, tasks, gamma)
    P_wt = wt_info["P"]

    L, K = genome.shape
    T = tasks.shape[1]
    dF = np.zeros((L * K, T), dtype=float)

    idx = 0
    for i in range(L):
        for j in range(K):
            mut = genome.copy()
            mut[i, j] = 1.0 - mut[i, j]
            P_mut = compute_performance(mut, tasks, gamma)["P"]
            dF[idx] = P_mut - P_wt
            idx += 1

    return dF


def mean_pairwise_tradeoff_modularity(
    genome: np.ndarray,
    tasks: np.ndarray,
    gamma: float,
) -> Tuple[float, float]:
    dF = mutation_effect_matrix(genome, tasks, gamma)
    T = dF.shape[1]

    tradeoffs = []
    modularities = []

    for a, b in combinations(range(T), 2):
        dFa = dF[:, a]
        dFb = dF[:, b]

        chi_ab = pairwise_tradeoff_for_two_tasks(dFa, dFb)
        M_ab = pairwise_modularity_for_two_tasks(dFa, dFb)

        if np.isfinite(chi_ab):
            tradeoffs.append(chi_ab)
        if np.isfinite(M_ab):
            modularities.append(M_ab)

    mean_tradeoff = np.nan if len(tradeoffs) == 0 else float(np.mean(tradeoffs))
    mean_modularity = np.nan if len(modularities) == 0 else float(np.mean(modularities))

    return mean_modularity, mean_tradeoff


# ============================================================
# 6. SNAPSHOT TRAJECTORY EXTRACTION
# ============================================================

def extract_snapshot_points_for_rep(
    rep: Dict,
    tasks: np.ndarray,
    gamma: float,
    n_target_points: int = 12,
) -> Tuple[np.ndarray, List[int]]:
    """
    Return
      pts: (n_target_points, 2)
      snap_steps: original saved snapshot step indices for this replicate

    If the replicate stalled early and fewer snapshots were saved,
    carry the last available point forward to n_target_points.
    """
    snap_steps = sorted(rep["snapshots"].keys())
    pts = []

    for step in snap_steps:
        genome = rep["snapshots"][step]["genome"]
        x, y = mean_pairwise_tradeoff_modularity(genome, tasks, gamma)
        pts.append([x, y])

    pts = np.array(pts, dtype=float)

    if len(pts) == 0:
        return np.full((n_target_points, 2), np.nan), snap_steps

    if len(pts) < n_target_points:
        last = pts[-1][None, :]
        pad = np.repeat(last, n_target_points - len(pts), axis=0)
        pts = np.vstack([pts, pad])
    elif len(pts) > n_target_points:
        pts = pts[:n_target_points]

    return pts, snap_steps


def extract_condition_trajectories(
    reps: List[Dict],
    tasks: np.ndarray,
    gamma: float,
    n_target_points: int = 12,
) -> Tuple[np.ndarray, List[List[int]]]:
    """
    Return
      trajs: (n_reps, n_target_points, 2)
      snapshot_steps_per_rep: original ragged step lists
    """
    out = []
    step_lists = []

    for rep in reps:
        traj, steps = extract_snapshot_points_for_rep(
            rep,
            tasks,
            gamma,
            n_target_points=n_target_points,
        )
        out.append(traj)
        step_lists.append(steps)

    return np.array(out, dtype=float), step_lists


# ============================================================
# 7. CONDITION-LEVEL GETTER WITH CACHE
# ============================================================

def get_or_compute_plane_trajectories(
    data_cfg: DataConfig,
    T: int,
    dT: float,
    m: int,
    alpha: float,
    force_recompute: bool = False,
) -> Optional[np.ndarray]:
    """
    Load plane trajectories from thin cache if present.
    Otherwise compute from saved snapshot genomes, then cache.

    Returns trajs with shape (n_reps, n_points, 2), or None if the
    underlying simulation cache is missing.
    """
    base_path = plane_cache_path(
        plane_cache_dir=data_cfg.plane_cache_dir,
        L=data_cfg.L,
        K=data_cfg.K,
        gamma=data_cfg.gamma,
        fitness_r=data_cfg.fitness_r,
        density=data_cfg.density,
        T=T,
        dT=dT,
        m=m,
        alpha=alpha,
    )

    if (not force_recompute) and os.path.exists(base_path + ".npz") and os.path.exists(base_path + "_meta.json"):
        trajs, _, _ = load_plane_condition(base_path)
        return trajs

    reps, _ = load_condition_if_exists(
        cache_dir=data_cfg.cache_dir,
        L=data_cfg.L,
        K=data_cfg.K,
        gamma=data_cfg.gamma,
        fitness_r=data_cfg.fitness_r,
        density=data_cfg.density,
        T=T,
        dT=dT,
        m=m,
        alpha=alpha,
    )
    if reps is None:
        return None

    tasks = load_task_map(
        task_cache_path(data_cfg.cache_dir, data_cfg.L, data_cfg.K,
                        data_cfg.gamma, data_cfg.fitness_r, T)
    )[dT]

    trajs, snapshot_steps_per_rep = extract_condition_trajectories(
        reps,
        tasks,
        data_cfg.gamma,
        n_target_points=12,
    )

    meta = {
        "L": data_cfg.L,
        "K": data_cfg.K,
        "gamma": data_cfg.gamma,
        "density": data_cfg.density,
        "T": T,
        "dT": dT,
        "m": m,
        "alpha": alpha,
        "n_reps": int(trajs.shape[0]),
        "n_points": int(trajs.shape[1]),
        "metric_definition": {
            "x": "mean pairwise Tikhonov-style modularity M over task pairs",
            "y": "mean pairwise Tikhonov-style tradeoff chi over task pairs",
        },
    }

    save_plane_condition(base_path, trajs, snapshot_steps_per_rep, meta)
    print(f"Cached plane coords: {base_path}.npz")

    return trajs


# ============================================================
# 8. PLOTTING
# ============================================================

def draw_arrow_on_last_segment(ax, xs, ys, color, lw, n_smooth: int = 4):
    """
    Draw an arrowhead at the end of the trajectory whose direction is
    estimated from the mean of up to `n_smooth` trailing displacement
    vectors. Carry-forward duplicates (zero displacement) are skipped so
    that stalled trajectories in row 2 still get a visible arrow.
    """
    if len(xs) < 2:
        return

    # collect finite, non-zero trailing displacements
    dx_list, dy_list = [], []
    x1, y1 = None, None
    for i in range(len(xs) - 2, -1, -1):
        x0, y0 = xs[i], ys[i]
        xi1, yi1 = xs[i + 1], ys[i + 1]
        if not np.isfinite([x0, y0, xi1, yi1]).all():
            continue
        ddx, ddy = xi1 - x0, yi1 - y0
        if abs(ddx) + abs(ddy) < 1e-12:
            continue  # skip carry-forward duplicates
        if x1 is None:
            x1, y1 = xi1, yi1  # arrowhead tip = last finite moving point
        dx_list.append(ddx)
        dy_list.append(ddy)
        if len(dx_list) >= n_smooth:
            break

    if x1 is None or len(dx_list) == 0:
        return

    # average direction over collected segments
    dx = float(np.mean(dx_list))
    dy = float(np.mean(dy_list))
    norm = np.sqrt(dx ** 2 + dy ** 2)
    if norm < 1e-12:
        return

    # phantom tail one step back along smoothed direction
    x0_phantom = x1 - dx / norm * norm
    y0_phantom = y1 - dy / norm * norm

    ax.annotate(
        "",
        xy=(x1, y1),
        xytext=(x0_phantom, y0_phantom),
        arrowprops=dict(
            arrowstyle="-|>",
            color=color,
            lw=lw * 0.8,
            mutation_scale=6,
            shrinkA=0.0,
            shrinkB=0.0,
        ),
        zorder=5,
    )


def plot_panel(
    ax,
    T: int,
    m: int,
    data_cfg: DataConfig,
    fig_cfg: FigureConfig,
    dt_color_map: Dict[float, Tuple[float, float, float]],
    force_recompute: bool = False,
):
    alpha_map = load_alpha_map(
        data_cfg.cache_dir, data_cfg.L, data_cfg.K, data_cfg.gamma,
        data_cfg.fitness_r, T
    )

    for dT in data_cfg.task_divs:
        if dT not in alpha_map:
            continue

        alpha = alpha_map[dT]

        trajs = get_or_compute_plane_trajectories(
            data_cfg=data_cfg,
            T=T,
            dT=dT,
            m=m,
            alpha=alpha,
            force_recompute=force_recompute,
        )
        if trajs is None:
            continue

        # faint replicate points, tinted with dT color
        color = dt_color_map[dT]
        for r in range(trajs.shape[0]):
            xs = trajs[r, :, 0]
            ys = trajs[r, :, 1]
            good = np.isfinite(xs) & np.isfinite(ys)
            if np.sum(good) < 1:
                continue

            ax.scatter(
                xs[good], ys[good],
                color=color,
                alpha=fig_cfg.faint_alpha,
                s=fig_cfg.faint_ms ** 2,
                linewidths=0,
                zorder=1,
            )

        # mean trajectory: line + arrow only (no markers)
        mean_x = np.nanmean(trajs[:, :, 0], axis=0)
        mean_y = np.nanmean(trajs[:, :, 1], axis=0)
        good = np.isfinite(mean_x) & np.isfinite(mean_y)
        if np.sum(good) < 2:
            continue

        ax.plot(
            mean_x[good], mean_y[good], "-",
            color=color,
            lw=fig_cfg.mean_lw,
            zorder=4,
        )


def make_figure(
    data_cfg: DataConfig,
    fig_cfg: FigureConfig,
    save_path: Optional[str] = None,
    force_recompute: bool = False,
):
    n_rows = 2
    n_cols = len(data_cfg.T_values)

    fig_w = fig_cfg.figsize_per_col[0] * n_cols
    fig_h = fig_cfg.figsize_per_col[1] * n_rows

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h))
    if n_cols == 1:
        axes = axes[:, np.newaxis]

    fig.subplots_adjust(
        left=0.09, right=0.97, top=0.92, bottom=0.10,
        wspace=0.18, hspace=0.18,
    )

    dt_color_map = make_dt_viridis_map(data_cfg.task_divs)
    panel_labels = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    pi = 0

    for col, T in enumerate(data_cfg.T_values):
        for row in range(n_rows):
            ax = axes[row, col]
            ax.set_box_aspect(1)

            if pi < len(panel_labels):
                ax.text(
                    -0.14, 1.08, panel_labels[pi],
                    transform=ax.transAxes,
                    fontsize=fig_cfg.panel_label_fs,
                    fontweight="bold",
                    ha="left", va="top",
                )
            pi += 1

            m = 1 if row == 0 else T

            plot_panel(
                ax=ax,
                T=T,
                m=m,
                data_cfg=data_cfg,
                fig_cfg=fig_cfg,
                dt_color_map=dt_color_map,
                force_recompute=force_recompute,
            )

            if row == 0:
                ax.set_title(f"T/K = {T / data_cfg.K:.1f}", fontsize=fig_cfg.title_fs, pad=6)

            if col == 0:
                ax.set_ylabel(
                    "Mutational tradeoff\n(mean pairwise $\\chi$)" +
                    ("\n$m = 1$" if row == 0 else "\n$m = T$")
                )
            else:
                ax.tick_params(labelleft=False)

            if row == n_rows - 1:
                ax.set_xlabel("Modularity $M$ (mean pairwise)")
            else:
                ax.tick_params(labelbottom=False)

            ax.set_xlim(0.0, 1.0)
            ax.set_ylim(-1.0, 1.0)
            ax.set_xticks([0.0, 0.5, 1.0])
            ax.set_yticks([-1.0, -0.5, 0.0, 0.5, 1.0])

    add_viridis_colorbar(fig, data_cfg.task_divs)

    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
        print(f"Saved: {save_path}")

    return fig


# ============================================================
# 9. CLI
# ============================================================

def parse_args():
    defaults = DataConfig()
    p = argparse.ArgumentParser()

    p.add_argument("--cache_dir", type=str, default=defaults.cache_dir)
    p.add_argument("--plane_cache_dir", type=str, default=defaults.plane_cache_dir)
    p.add_argument("--save_dir", type=str, default="figures_out")
    p.add_argument("--filename", type=str, default="figS_tradeoff_modularity")
    p.add_argument("--fmt", type=str, default="pdf")

    p.add_argument("--L", type=int, default=defaults.L)
    p.add_argument("--K", type=int, default=defaults.K)
    p.add_argument("--gamma", type=float, default=defaults.gamma)
    p.add_argument("--fitness_r", type=float, default=defaults.fitness_r)
    p.add_argument("--density", type=float, default=defaults.density)

    p.add_argument("--T", type=int, nargs="+", default=defaults.T_values, dest="T_values")
    p.add_argument("--dT", type=float, nargs="+", default=defaults.task_divs, dest="task_divs")

    p.add_argument("--force_recompute", action="store_true")
    p.add_argument("--no_show", action="store_true")
    return p.parse_args()


# ============================================================
# 10. MAIN
# ============================================================

if __name__ == "__main__":
    args = parse_args()

    data_cfg = DataConfig(
        cache_dir=args.cache_dir,
        plane_cache_dir=args.plane_cache_dir,
        L=args.L,
        K=args.K,
        gamma=args.gamma,
        fitness_r=args.fitness_r,
        density=args.density,
        T_values=args.T_values,
        task_divs=args.task_divs,
    )

    fig_cfg = FigureConfig()

    os.makedirs(args.save_dir, exist_ok=True)
    os.makedirs(data_cfg.plane_cache_dir, exist_ok=True)

    save_path = os.path.join(
        args.save_dir,
        f"{args.filename}"
        f"_gamma{data_cfg.gamma}_fr{data_cfg.fitness_r}"
        f"_density{data_cfg.density:.4f}.{args.fmt}"
    )

    fig = make_figure(
        data_cfg=data_cfg,
        fig_cfg=fig_cfg,
        save_path=save_path,
        force_recompute=args.force_recompute,
    )

    if args.no_show:
        plt.close(fig)
    else:
        plt.show()