#!/usr/bin/env python3
"""
compare_K.py
============
Differentiation and optimization against the task-to-program ratio T/K, with
one line per program number.

The argument that program number sets the limit under simultaneous selection
requires comparing caches at matched T/K, not matched T. At K=4 the grid
T = 2, 4, 6, 8 covers T/K = 0.5, 1, 1.5, 2; the same T values at K=6 cover only
0.33 to 1.33, so a K=4 versus K=6 comparison at equal T conflates "more tasks
than programs" with "more tasks". Plotting against T/K makes the mismatch
visible and shows how much of the informative range each cache actually spans.

  Rows     differentiation, optimization
  Columns  sequential (m=1), simultaneous (m=T)
  x-axis   T/K, with the boundary at 1 marked
  Lines    one per task divergence (colour), one style per K

Under simultaneous selection the two K curves should coincide when plotted this
way if program number is what sets the limit. Under sequential selection the
decline should begin left of T/K = 1 for both.

Usage:
  python3 compare_K.py --K 4 6
  python3 compare_K.py --K 4 6 --dT 0.2 1.4 --cutoff 200
  python3 compare_K.py --K 4 6 8 --no_plot
"""

import argparse
import os
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np

import figlib as FL


LINESTYLES = ['-', '--', ':', '-.']


def collect(K_values: List[int], T_by_K: Dict[int, List[int]], args,
            cutoff: FL.Cutoff):
    """{K: {'spec': CacheSpec, 'data': grid}} for every K with a cache."""
    out = {}
    for K in K_values:
        spec = FL.CacheSpec(cache_dir=args.cache_dir, L=args.L, K=K,
                            gamma=args.gamma, fitness_r=args.fitness_r,
                            density=args.density,
                            T_values=T_by_K.get(K, args.T_values),
                            task_divs=args.task_divs)
        print(f'Loading K={K} ...')
        data = FL.load_grid(spec, verbose=False)
        found = {T: sorted(data.get(T, {}).get(spec.task_divs[0], {}))
                 for T in spec.T_values}
        found = {T: ms for T, ms in found.items() if ms}
        if not found:
            print(f'  no conditions found for K={K}')
            continue
        print(f'  T present: {sorted(found)}')
        out[K] = {'spec': spec, 'data': data}
    return out


def series(entry, metric: str, dT: float, regime: str, cutoff: FL.Cutoff):
    """(T/K, mean, sd) across the T values present for one K."""
    spec, data = entry['spec'], entry['data']
    xs, ys, sds = [], [], []
    for T in spec.T_values:
        by_m = data.get(T, {}).get(dT, {})
        m = FL.resolve_m('min' if regime == 'sequential' else 'T',
                         by_m.keys(), T)
        if m is None:
            continue
        _, vals = FL.metric_values(by_m[m], metric, cutoff)
        mu, sd = FL.mean_sd(vals)
        if np.isfinite(mu):
            xs.append(T / spec.K); ys.append(mu); sds.append(sd)
    return np.array(xs), np.array(ys), np.array(sds)


def print_table(caches, cutoff, task_divs):
    print(f'\n{"=" * 92}')
    print(f'DIFFERENTIATION AND OPTIMIZATION BY T/K  ({cutoff.label()})')
    print('Rows are matched on T/K, so a K=4 row at T/K=2 (T=8) is comparable '
          'to a K=6 row at T/K=2 (T=12).')
    print('=' * 92)
    head = (f'{"K":>3} {"T":>3} {"T/K":>5} {"dT":>5} '
            f'{"diff_m1":>16} {"diff_mT":>16} {"opt_m1":>16} {"opt_mT":>16}')
    print(head)
    print('-' * len(head))

    for K in sorted(caches):
        spec, data = caches[K]['spec'], caches[K]['data']
        for T in spec.T_values:
            for dT in task_divs:
                by_m = data.get(T, {}).get(dT, {})
                if not by_m:
                    continue
                cells = []
                for metric in ('differentiation', 'optimization'):
                    for sel in ('min', 'T'):
                        m = FL.resolve_m(sel, by_m.keys(), T)
                        if m is None:
                            cells.append('--')
                            continue
                        _, v = FL.metric_values(by_m[m], metric, cutoff)
                        mu, sd = FL.mean_sd(v)
                        cells.append(f'{mu:.4f}+/-{sd:.4f}')
                print(f'{K:>3} {T:>3} {T / K:>5.2f} {dT:>5.1f} '
                      f'{cells[0]:>16} {cells[1]:>16} '
                      f'{cells[2]:>16} {cells[3]:>16}')


def make_figure(caches, task_divs, cutoff, save_path=None):
    FL.apply_style()
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 8.0), squeeze=False)
    fig.subplots_adjust(hspace=0.30, wspace=0.22,
                        left=0.12, right=0.88, top=0.90, bottom=0.14)

    colors = FL.dt_colors(task_divs)
    styles = {K: LINESTYLES[i % len(LINESTYLES)]
              for i, K in enumerate(sorted(caches))}
    rows = [('differentiation', 'Degree of differentiation'),
            ('optimization', 'Degree of optimization')]
    cols = [('sequential', 'Sequential ($m=1$)'),
            ('simultaneous', 'Simultaneous ($m=T$)')]

    for r, (metric, ylabel) in enumerate(rows):
        for c, (regime, title) in enumerate(cols):
            ax = axes[r][c]
            ax.set_box_aspect(1)
            ax.text(-0.15, 1.08, FL.panel_label(r * 2 + c),
                    transform=ax.transAxes, fontsize=14, fontweight='bold',
                    va='top', ha='left')
            if r == 0:
                ax.set_title(title, fontsize=12, pad=6)

            for K in sorted(caches):
                for dT in task_divs:
                    xs, ys, sds = series(caches[K], metric, dT, regime, cutoff)
                    if xs.size == 0:
                        continue
                    ax.errorbar(xs, ys, yerr=sds, fmt='o', ls=styles[K],
                                color=colors[dT], lw=0.75, ms=4,
                                markerfacecolor='none', markeredgecolor=colors[dT],
                                capsize=2, capthick=0.8, elinewidth=0.8)

            ax.axvline(1.0, color='gray', ls=':', lw=1.0, alpha=0.8, zorder=0)
            if r == 0:
                ax.annotate('$T=K$', xy=(1.0, 1.0),
                            xycoords=('data', 'axes fraction'),
                            xytext=(3, -3), textcoords='offset points',
                            fontsize=8, color='gray', ha='left', va='top')
            ax.axhline(1, color='gray', ls='--', lw=0.8, alpha=0.5)
            ax.set_ylim(0, 1.05)
            ax.set_xlabel('Tasks per program, $T/K$')
            if c == 0:
                ax.set_ylabel(ylabel)
            else:
                ax.tick_params(labelleft=False)

    handles = [plt.Line2D([], [], color='0.3', ls=styles[K], label=f'$K={K}$')
               for K in sorted(caches)]
    axes[0][1].legend(handles=handles, fontsize=9, frameon=False, loc='best')

    FL.add_dt_colorbar(fig, task_divs)
    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
        print(f'\nSaved: {save_path}')
    return fig


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--cache_dir', default='simulation_cache')
    p.add_argument('--save_dir', default='figures_out')
    p.add_argument('--filename', default='compare_K')
    p.add_argument('--fmt', default='pdf')
    p.add_argument('--L', type=int, default=100)
    p.add_argument('--K', type=int, nargs='+', default=[4, 6], dest='K_values')
    p.add_argument('--T', type=int, nargs='+', default=[2, 3, 4, 6, 8, 9, 12, 16],
                   dest='T_values',
                   help='Candidate T values; those without a cache are skipped.')
    p.add_argument('--dT', type=float, nargs='+', default=[0.2, 0.8, 1.4],
                   dest='task_divs')
    p.add_argument('--gamma', type=float, default=1.0)
    p.add_argument('--fitness_r', type=float, default=0.0)
    p.add_argument('--density', type=float, default=0.25)
    p.add_argument('--cutoff', type=int, default=200)
    p.add_argument('--cutoff_kind', default='substitutions',
                   choices=['substitutions', 'exposure'])
    p.add_argument('--no_plot', action='store_true')
    p.add_argument('--no_show', action='store_true')
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    cutoff = FL.Cutoff(args.cutoff_kind, args.cutoff)

    caches = collect(args.K_values, {}, args, cutoff)
    if not caches:
        raise SystemExit('No caches found.')

    print_table(caches, cutoff, args.task_divs)

    if not args.no_plot:
        os.makedirs(args.save_dir, exist_ok=True)
        path = os.path.join(
            args.save_dir,
            f'{args.filename}_sub{args.cutoff}_gamma{args.gamma}'
            f'_fr{args.fitness_r}_density{args.density:.4f}.{args.fmt}')
        fig = make_figure(caches, args.task_divs, cutoff, save_path=path)
        if args.no_show:
            plt.close(fig)
        else:
            plt.show()