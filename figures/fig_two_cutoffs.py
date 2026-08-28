#!/usr/bin/env python3
"""
fig_two_cutoffs.py
==================
One condition shown at two points along the substitution trajectory, in a
single figure.

  Rows     differentiation, optimization
  Columns  sequential (m=1), simultaneous (m=T)
  x-axis   number of tasks, with T=K marked
  Lines    one colour per task divergence; solid for the earlier comparison
           point, dashed for the later one

This is for analyses in which a parameter change alters how quickly populations
reach an absorbing state rather than the state they reach. Comparing at one
substitution count then catches the affected populations part-way, and the
difference between the solid and dashed lines is the size of that shortfall.
Where the dashed lines recover the pattern seen in the main figures, the effect
of the parameter is on rate rather than on outcome.

Only a few task divergences are shown, because two comparison points across the
full grid produces too many lines to read.

Usage:
  python3 fig_two_cutoffs.py --density 0.5 --filename FS5_density
  python3 fig_two_cutoffs.py --K 6 --cutoffs 200 400 --filename FS_K6_cutoffs
  python3 fig_two_cutoffs.py --density 0.5 --dT 0.2 0.8 1.4
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
from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as np

import figlib as FL


REGIMES = [('min', 'Sequential ($m=1$)'),
           ('T', 'Simultaneous ($m=T$)')]

METRICS = [('differentiation', 'Degree of differentiation'),
           ('optimization', 'Degree of optimization')]

LINESTYLES = ['-', '--', ':', '-.']


def make_figure(data, spec: FL.CacheSpec, cutoffs: List[FL.Cutoff],
                show_K_line: bool = True, save_path: Optional[str] = None):
    FL.apply_style()
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 8.0), squeeze=False)
    fig.subplots_adjust(hspace=0.30, wspace=0.22,
                        left=0.12, right=0.94, top=0.90, bottom=0.16)

    colors = FL.dt_colors(spec.task_divs)
    t_values = np.array(spec.T_values)

    for r, (metric, ylabel) in enumerate(METRICS):
        for c, (selector, title) in enumerate(REGIMES):
            ax = axes[r][c]
            ax.set_box_aspect(1)
            ax.text(-0.15, 1.08, FL.panel_label(r * 2 + c),
                    transform=ax.transAxes, fontsize=14, fontweight='bold',
                    va='top', ha='left')
            if r == 0:
                ax.set_title(title, fontsize=12, pad=6)

            for k, cutoff in enumerate(cutoffs):
                for dT in spec.task_divs:
                    xs, ys, sds = [], [], []
                    for T in spec.T_values:
                        by_m = data.get(T, {}).get(dT, {})
                        m = FL.resolve_m(selector, by_m.keys(), T)
                        if m is None:
                            continue
                        _, vals = FL.metric_values(by_m[m], metric, cutoff)
                        mu, sd = FL.mean_sd(vals)
                        if np.isfinite(mu):
                            xs.append(T); ys.append(mu); sds.append(sd)
                    if not xs:
                        continue
                    ax.errorbar(xs, ys, yerr=sds, fmt='o',
                                ls=LINESTYLES[k % len(LINESTYLES)],
                                color=colors[dT], lw=0.9, ms=4,
                                markerfacecolor='none',
                                markeredgecolor=colors[dT],
                                capsize=2, capthick=0.8, elinewidth=0.8)

            ax.set_xticks(t_values)
            ax.set_xticklabels([str(int(v)) for v in t_values])
            ax.set_xlabel('Number of tasks')

            if show_K_line and t_values.min() <= spec.K <= t_values.max():
                ax.axvline(spec.K, color='gray', ls=':', lw=1.0, alpha=0.8,
                           zorder=0)
                if r == 0:
                    ax.annotate(f'$T=K={spec.K}$', xy=(spec.K, 1.0),
                                xycoords=('data', 'axes fraction'),
                                xytext=(3, -3), textcoords='offset points',
                                fontsize=8, color='gray', ha='left', va='top')

            ax.axhline(1, color='gray', ls='--', lw=0.8, alpha=0.5)
            ax.set_ylim(0, 1.05)
            if c == 0:
                ax.set_ylabel(ylabel)
            else:
                ax.tick_params(labelleft=False)

    handles = [plt.Line2D([], [], color='0.3',
                          ls=LINESTYLES[k % len(LINESTYLES)], marker='o',
                          markerfacecolor='none',
                          label=f'{cutoff.value} substitutions')
               for k, cutoff in enumerate(cutoffs)]
    axes[0][1].legend(handles=handles, fontsize=9, frameon=False, loc='best')

    FL.add_dt_colorbar(fig, spec.task_divs)
    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
        print(f'Saved: {save_path}')
    return fig


def print_summary(data, spec: FL.CacheSpec, cutoffs):
    print(f'\n{"=" * 104}')
    print(f'SUMMARY  {spec.label()}')
    print('Each condition at both comparison points.')
    print('=' * 104)
    head = (f'{"T":>3} {"dT":>5} {"cut":>5} '
            f'{"diff_m1":>16} {"diff_mT":>16} {"opt_m1":>16} {"opt_mT":>16}')
    print(head)
    print('-' * len(head))

    for T in spec.T_values:
        for dT in spec.task_divs:
            by_m = data.get(T, {}).get(dT, {})
            if not by_m:
                continue
            for cutoff in cutoffs:
                cells = []
                for metric, _ in METRICS:
                    for selector, _ in REGIMES:
                        m = FL.resolve_m(selector, by_m.keys(), T)
                        if m is None:
                            cells.append('--')
                            continue
                        _, v = FL.metric_values(by_m[m], metric, cutoff)
                        mu, sd = FL.mean_sd(v)
                        cells.append(f'{mu:.4f}+/-{sd:.4f}')
                print(f'{T:>3} {dT:>5.1f} {cutoff.value:>5} '
                      f'{cells[0]:>16} {cells[1]:>16} '
                      f'{cells[2]:>16} {cells[3]:>16}')


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--cache_dir', default=_repo_path('simulation_cache'))
    p.add_argument('--save_dir', default=_repo_path('figures_out'))
    p.add_argument('--filename', default='FS_two_cutoffs')
    p.add_argument('--fmt', default='pdf')
    p.add_argument('--L', type=int, default=100)
    p.add_argument('--K', type=int, default=4)
    p.add_argument('--gamma', type=float, default=1.0)
    p.add_argument('--fitness_r', type=float, default=0.0)
    p.add_argument('--density', type=float, default=0.5)
    p.add_argument('--T', type=int, nargs='+', default=[2, 4, 6, 8],
                   dest='T_values')
    p.add_argument('--dT', type=float, nargs='+', default=[0.2, 0.8, 1.4],
                   dest='task_divs',
                   help='Few values only; two comparison points across the '
                        'full grid is unreadable.')
    p.add_argument('--cutoffs', type=int, nargs='+', default=[200, 400])
    p.add_argument('--no_K_line', action='store_true')
    p.add_argument('--no_show', action='store_true')
    p.add_argument('--no_summary', action='store_true')
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()

    spec = FL.CacheSpec(cache_dir=args.cache_dir, L=args.L, K=args.K,
                        gamma=args.gamma, fitness_r=args.fitness_r,
                        density=args.density, T_values=args.T_values,
                        task_divs=args.task_divs)
    cutoffs = [FL.Cutoff('substitutions', c) for c in args.cutoffs]

    print(f'Loading {spec.label()} ...')
    data = FL.load_grid(spec)

    os.makedirs(args.save_dir, exist_ok=True)
    tag = '_'.join(str(c) for c in args.cutoffs)
    path = os.path.join(
        args.save_dir,
        f'{args.filename}_cut{tag}_gamma{args.gamma}_fr{args.fitness_r}'
        f'_K{args.K}_density{args.density:.4f}.{args.fmt}')

    fig = make_figure(data, spec, cutoffs,
                      show_K_line=not args.no_K_line, save_path=path)
    if not args.no_summary:
        print_summary(data, spec, cutoffs)
    if args.no_show:
        plt.close(fig)
    else:
        plt.show()