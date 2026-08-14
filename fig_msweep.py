#!/usr/bin/env python3
"""
fig_msweep.py
=============
Gain over sequential selection as a function of simultaneity m, for both
response metrics.

  Rows    differentiation and optimization
  Columns one per task count T
  x-axis  m, the number of task-specific phenotypes contributing jointly
          to fitness at each fixation event
  y-axis  gain(m) = metric(m) - metric(m=1)
  Lines   one per task divergence, viridis_r

--------------------------------------------------------------------
Pairing
--------------------------------------------------------------------
Replicate i uses the same initial genome and the same task ensemble at every m
(see simulate.py, "Replicate structure"). The contrast between m and m=1 is
therefore PAIRED, and gain(m) is computed as the mean of within-replicate
differences with the standard error of those differences,

    SE = s_d / sqrt(n),    d_i = X_i(m) - X_i(1),

rather than as a difference of two independent means. Pairing removes the
ensemble and initial-genome contributions from the error and requires no
independence assumption between conditions.

m=1 is the baseline by construction, so its gain and error are exactly zero.

--------------------------------------------------------------------
Framing
--------------------------------------------------------------------
gain(m) measures how much differentiation or optimization a lineage with
historical simultaneity m has built up relative to purely sequential selection.
It is anchored at zero at the most sequential case and grows with m, which is
the orientation needed to ask whether intermediate m acts as a precursor to a
transition to full simultaneity.

The complementary quantity, the deficit relative to m=T, is obtained by setting
baseline='T'; the two describe the same data anchored at opposite ends, since
gain(m) + deficit(m) is constant for fixed T and dT.

`fraction_of_gain` reports 100 * gain(m) / gain(T), the percentage of the fully
simultaneous gain already present at intermediate m. It is undefined when
gain(T) is not clearly positive, and is reported as nan rather than as a large
or negative percentage in that case.

Usage:
  python3 fig_msweep.py
  python3 fig_msweep.py --T 4 8 --cutoffs 200
  python3 fig_msweep.py --cutoff_kind exposure --cutoffs 50
"""

import argparse
import os
from dataclasses import dataclass, field
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

def make_figure(data, spec: FL.CacheSpec, cutoff: FL.Cutoff,
                fig_cfg: FigConfig, save_path: Optional[str] = None):
    FL.apply_style()
    n_cols = len(spec.T_values)
    fig, axes = plt.subplots(
        2, n_cols, squeeze=False,
        figsize=(fig_cfg.figsize_per_col[0] * n_cols,
                 fig_cfg.figsize_per_col[1] * 2))
    fig.subplots_adjust(hspace=0.28, wspace=0.15,
                        left=0.10, right=0.94, top=0.92, bottom=0.10)

    colors = FL.dt_colors(spec.task_divs)
    rows = [('differentiation', r'$\Delta$ degree of differentiation'),
            ('optimization', r'$\Delta$ degree of optimization')]

    for col, T in enumerate(spec.T_values):
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
                ax.errorbar(ms, gains, yerr=ses, fmt='-o', color=colors[dT],
                            lw=fig_cfg.line_width, ms=fig_cfg.marker_size,
                            markerfacecolor='none', markeredgecolor=colors[dT],
                            capsize=2, capthick=1.0, elinewidth=1.0, zorder=3)

            ax.axhline(0, color='gray', ls='--', lw=0.8, alpha=0.5)
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

    FL.add_dt_colorbar(fig, spec.task_divs, orientation='vertical')
    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
        print(f'Saved: {save_path}')
    return fig


# ============================================================
# SUMMARY
# ============================================================

def print_summary(data, spec: FL.CacheSpec, cutoff: FL.Cutoff,
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
    p.add_argument('--cache_dir', default='simulation_cache')
    p.add_argument('--save_dir', default='figures_out')
    p.add_argument('--filename', default='F4')
    p.add_argument('--fmt', default='pdf')
    p.add_argument('--L', type=int, default=100)
    p.add_argument('--K', type=int, default=4)
    p.add_argument('--gamma', type=float, default=1.0)
    p.add_argument('--fitness_r', type=float, default=0.0)
    p.add_argument('--density', type=float, default=0.25)
    p.add_argument('--T', type=int, nargs='+', default=[2, 4, 8],
                   dest='T_values')
    p.add_argument('--dT', type=float, nargs='+',
                   default=[0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4],
                   dest='task_divs')
    p.add_argument('--cutoff', type=int, default=200)
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
    data = FL.load_grid(spec)

    os.makedirs(args.save_dir, exist_ok=True)
    tag = ('' if args.cutoff_kind == 'substitutions'
           else f'_exposure{args.exposure_mode}')
    path = os.path.join(
        args.save_dir,
        f'{args.filename}{tag}_sub{args.cutoff}'
        f'_gamma{args.gamma}_fr{args.fitness_r}'
        f'_K{args.K}_density{args.density:.4f}.{args.fmt}')

    fig = make_figure(data, spec, cutoff, fig_cfg, save_path=path)
    if not args.no_summary:
        print_summary(data, spec, cutoff, fig_cfg)
    if args.no_show:
        plt.close(fig)
    else:
        plt.show()
