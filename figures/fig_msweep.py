#!/usr/bin/env python3
"""
fig_msweep.py
=============
Figure S6: gain over sequential selection, metric(m) - metric(m=1), against m.

Bands are the standard error of the within-replicate paired difference, not the
spread across replicates. Replicate i shares its initial genome and task
ensemble across every m, so the paired difference removes the between-world
variance that dominates the error bands of Figure 4; the bands here are
correspondingly much narrower, and answer a different question.

Colour is task divergence, on the scale shared with every other figure, so this
figure carries no colour key of its own.

One task number is loaded at a time. This is the only figure that needs every
simultaneity level, and the whole grid does not fit in memory at once.
"""

# --- repo root on sys.path, so this script runs from any working directory ---
import os as _os
import sys as _sys

_REPO_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
if _REPO_ROOT not in _sys.path:
    _sys.path.insert(0, _REPO_ROOT)


def _repo_path(*parts):
    """Path anchored at the repository root, independent of the caller's cwd."""
    return _os.path.join(_REPO_ROOT, *parts)

# ---------------------------------------------------------------------------


import argparse
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

import figlib as FL


MIN_DENOMINATOR = 1e-3      # below this, fraction-of-gain is meaningless


@dataclass
class FigConfig:
    baseline: object = 1                # 1 = sequential; 'T' for the deficit view
    figsize_per_col: Tuple[float, float] = (4.5, 4.5)
    line_width: float = 0.75
    marker_size: float = 5.0
    title_fontsize: int = 14
    label_fontsize: int = 12


# ============================================================
# CORE
# ============================================================

def gain_curve(by_m: Dict[int, List], metric: str, cutoff: FL.Cutoff,
               baseline_m: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray,
                                         np.ndarray]:
    """(m, gain, se, n) swept over the available simultaneity levels.

    The baseline level itself is included with gain and se exactly zero; every
    other level is a paired contrast against it.
    """
    if baseline_m not in by_m:
        return (np.array([]),) * 4

    idx_b, vals_b = FL.metric_values(by_m[baseline_m], metric, cutoff)
    if idx_b.size == 0:
        return (np.array([]),) * 4

    ms, gains, ses, ns = [], [], [], []
    for m in sorted(by_m):
        if m == baseline_m:
            ms.append(m); gains.append(0.0); ses.append(0.0)
            ns.append(int(idx_b.size))
            continue
        idx_a, vals_a = FL.metric_values(by_m[m], metric, cutoff)
        g, se, n = FL.paired_difference(idx_a, vals_a, idx_b, vals_b)
        if np.isfinite(g):
            ms.append(m); gains.append(g); ses.append(se); ns.append(n)

    return (np.array(ms, dtype=float), np.array(gains), np.array(ses),
            np.array(ns, dtype=int))


def fraction_of_gain(ms: np.ndarray, gains: np.ndarray,
                     T: int) -> Dict[int, float]:
    """Percentage of the fully simultaneous gain present at each m.

    Returns nan where the m=T gain is too small for the ratio to carry meaning;
    a percentage of a near-zero denominator is noise, not a result.
    """
    if ms.size == 0 or T not in ms.astype(int):
        return {}
    denom = float(gains[np.flatnonzero(ms.astype(int) == T)[0]])
    if not np.isfinite(denom) or denom < MIN_DENOMINATOR:
        return {int(m): np.nan for m in ms}
    return {int(m): 100.0 * float(g) / denom for m, g in zip(ms, gains)}


# ============================================================
# FIGURE
# ============================================================

def make_figure(load_T, spec: FL.CacheSpec, cutoff: FL.Cutoff,
                fig_cfg: FigConfig, save_path: Optional[str] = None):
    """`load_T` is called once per task number and its result released before
    the next. This figure is the only one that needs every simultaneity level,
    and holding the whole grid at once does not fit in memory."""
    FL.apply_style()
    n_cols = len(spec.T_values)
    fig, axes = plt.subplots(
        2, n_cols, squeeze=False,
        figsize=(fig_cfg.figsize_per_col[0] * n_cols,
                 fig_cfg.figsize_per_col[1] * 2))
    fig.subplots_adjust(hspace=0.28, wspace=0.15,
                        left=0.10, right=0.94, top=0.92, bottom=0.10)

    colors = FL.dt_colors(spec.task_divs)
    rows = [('differentiation',
             'Gain in differentiation\nover sequential selection'),
            ('optimization',
             'Gain in optimization\nover sequential selection')]

    for col, T in enumerate(spec.T_values):
        data = load_T(T)
        for row, (metric, ylabel) in enumerate(rows):
            ax = axes[row][col]
            ax.set_box_aspect(1)
            if row == 0:
                ax.set_title(f'Number of tasks = {T}',
                             fontsize=fig_cfg.title_fontsize, pad=8)

            seen_m = set()
            for dT in spec.task_divs:
                by_m = data.get(T, {}).get(dT, {})
                if not by_m:
                    continue
                base = FL.resolve_m(fig_cfg.baseline, by_m.keys(), T)
                if base is None:
                    continue
                ms, gains, ses, _ = gain_curve(by_m, metric, cutoff, base)
                if ms.size == 0:
                    continue
                seen_m.update(ms.tolist())
                # Bands are +/- 1 SE of the within-replicate paired
                # difference here, not the SD across replicates.
                FL.band(ax, ms, gains, ses, color=colors[dT], ls=FL.LS_MT)

            ax.axhline(0, color='0.68', ls='-', lw=0.9, zorder=0)
            if seen_m:
                ax.set_xticks(sorted(seen_m))
            if row == 1:
                ax.set_xlabel('Number of tasks under selection ($m$)',
                              fontsize=fig_cfg.label_fontsize)
            if col == 0:
                ax.set_ylabel(ylabel, fontsize=fig_cfg.label_fontsize)
            else:
                ax.tick_params(labelleft=False)

    # shared y-limits per row so columns are directly comparable
    for row in range(2):
        lims = [axes[row][c].get_ylim() for c in range(n_cols)]
        lo, hi = min(l[0] for l in lims), max(l[1] for l in lims)
        for c in range(n_cols):
            axes[row][c].set_ylim(lo, hi)

    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
        print(f'Saved: {save_path}')
    return fig


# ============================================================
# SUMMARY
# ============================================================

def print_summary(load_T, spec: FL.CacheSpec, cutoff: FL.Cutoff,
                  fig_cfg: FigConfig):
    print(f'\n{"=" * 104}')
    print(f'SUMMARY  gain over m={fig_cfg.baseline}, paired  '
          f'({cutoff.label()})  {spec.label()}')
    print('SE is the standard error of within-replicate differences; n is the '
          'number of matched replicates.')
    print('%gain is the percentage of the m=T gain present at that m, nan '
          'where the m=T gain is too small to divide by.')
    print('=' * 104)
    head = (f'{"T":>3} {"dT":>5} {"m":>3} {"n":>5} '
            f'{"gain_diff":>10} {"SE":>8} {"%gain":>7}   '
            f'{"gain_opt":>10} {"SE":>8} {"%gain":>7}')
    print(head)
    print('-' * len(head))

    for T in spec.T_values:
        data = load_T(T)
        for dT in spec.task_divs:
            by_m = data.get(T, {}).get(dT, {})
            if not by_m:
                continue
            base = FL.resolve_m(fig_cfg.baseline, by_m.keys(), T)
            if base is None:
                continue

            md, gd, sd, nd = gain_curve(by_m, 'differentiation', cutoff, base)
            mo, go, so, no = gain_curve(by_m, 'optimization', cutoff, base)
            fd = fraction_of_gain(md, gd, T)
            fo = fraction_of_gain(mo, go, T)

            d_lookup = {int(m): (g, s, n) for m, g, s, n in zip(md, gd, sd, nd)}
            o_lookup = {int(m): (g, s, n) for m, g, s, n in zip(mo, go, so, no)}

            for m in sorted(set(d_lookup) | set(o_lookup)):
                dg, ds, dn = d_lookup.get(m, (np.nan, np.nan, 0))
                og, os_, on = o_lookup.get(m, (np.nan, np.nan, 0))
                print(f'{T:>3} {dT:>5.1f} {m:>3} {max(dn, on):>5} '
                      f'{dg:>10.4f} {ds:>8.4f} {fd.get(m, np.nan):>7.1f}   '
                      f'{og:>10.4f} {os_:>8.4f} {fo.get(m, np.nan):>7.1f}')


# ============================================================
# CLI
# ============================================================

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--cache_dir', default=_repo_path('simulation_cache'))
    p.add_argument('--save_dir', default=_repo_path('figures_out'))
    p.add_argument('--filename', default='F4')
    p.add_argument('--fmt', default='pdf')
    p.add_argument('--L', type=int, default=100)
    p.add_argument('--K', type=int, default=4)
    p.add_argument('--gamma', type=float, default=1.0)
    p.add_argument('--fitness_r', type=float, default=0.0)
    p.add_argument('--density', type=float, default=0.25)
    p.add_argument('--T', type=int, nargs='+', default=[2, 4, 6, 8],
                   dest='T_values')
    p.add_argument('--dT', type=float, nargs='+',
                   default=[0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4],
                   dest='task_divs')
    p.add_argument('--cutoff', type=int, default=400)
    p.add_argument('--cutoff_kind', default='substitutions',
                   choices=['substitutions', 'exposure'])
    p.add_argument('--exposure_mode', default='realized',
                   choices=['realized', 'expected'])
    p.add_argument('--no_show', action='store_true')
    p.add_argument('--no_summary', action='store_true')
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()

    spec = FL.CacheSpec(cache_dir=args.cache_dir, L=args.L, K=args.K,
                        gamma=args.gamma, fitness_r=args.fitness_r,
                        density=args.density, T_values=args.T_values,
                        task_divs=args.task_divs)
    cutoff = FL.Cutoff(args.cutoff_kind, args.cutoff, args.exposure_mode)
    fig_cfg = FigConfig()

    print(f'Loading {spec.label()} ...')
    from dataclasses import replace as _replace

    def load_T(T):
        return FL.load_grid(_replace(spec, T_values=[T]), verbose=False)

    os.makedirs(args.save_dir, exist_ok=True)
    tag = ('' if args.cutoff_kind == 'substitutions'
           else f'_exposure{args.exposure_mode}')
    path = os.path.join(
        args.save_dir,
        f'{args.filename}{tag}_sub{args.cutoff}'
        f'_gamma{args.gamma}_fr{args.fitness_r}'
        f'_K{args.K}_density{args.density:.4f}.{args.fmt}')

    fig = make_figure(load_T, spec, cutoff, fig_cfg, save_path=path)
    if not args.no_summary:
        print_summary(load_T, spec, cutoff, fig_cfg)
    if args.no_show:
        plt.close(fig)
    else:
        plt.show()