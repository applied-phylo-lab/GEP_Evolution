#!/usr/bin/env python3
"""
compare_K.py
============
Differentiation and optimization against the task-to-program ratio T/K, with
one line style per program number.

  Rows     differentiation, optimization
  Columns  sequential (m=1), simultaneous (m=T)
  x-axis   T/K, with the boundary T=K marked at 1
  Lines    one per task divergence (colour), one style per K

--------------------------------------------------------------------
Why T/K
--------------------------------------------------------------------
The claim under test is that program number sets where differentiation begins
to fail under simultaneous selection, and does not set it under sequential
selection. That is a claim about the POSITION of the decline, so the axis has
to be the one on which the position should be fixed: tasks per program.

Comparing program numbers at equal T instead would confound having more tasks
than programs with simply having more tasks, since T=8 is twice the boundary at
K=4 but only 1.33 times it at K=6.

The two curves are not expected to coincide numerically. L is held fixed while
K varies, so a larger K means a larger genotype matrix and more mutational
targets, and T task optima in L-dimensional space do not rescale geometrically
with T. What should coincide is the shape: a decline that sets in around T/K=1
under simultaneous selection for every K, and one that begins well left of it
under sequential selection.

Both program numbers must span the same range of T/K for the overlay to mean
anything. At K=4 the grid T = 2, 4, 6, 8 gives 0.5, 1, 1.5, 2; matching that at
K=6 requires T = 3, 6, 9, 12. `run_batch.py --specs programs --K 6 --T 3 6 9 12`
generates the missing points.

`--x_axis tasks` plots against T instead, with T=K marked separately for each
program number. That view answers a different question -- what adding programs
buys at a fixed task number -- which `print_interaction` reports numerically.

--------------------------------------------------------------------
Why the cutoff scales with K
--------------------------------------------------------------------
The genotype matrix has L*K entries, so a fixed substitution count is a
different fraction of the available change at each program number: 200
substitutions is half of L*K at K=4 but a quarter at K=8. Comparing program
numbers at a fixed substitution count therefore has the same defect as
comparing them at a fixed task number -- it holds the wrong quantity constant.

The default `--cutoff_scale genome` scales the cutoff with K, so that every
program number is evaluated after the same number of substitutions per genotype
entry: with a base of 200 at the smallest K in the comparison, K = 4, 6, 8 are
read at 200, 300 and 400 substitutions. `--cutoff_scale fixed` applies one
cutoff to every K.

In practice this matters only for the sequential panels. Under simultaneous
selection nearly every replicate reaches an absorbing state well before any of
these cutoffs, so those panels are effectively cutoff-independent; under
sequential selection the trajectories are still climbing and the choice is
material.

Usage:
  python3 compare_K.py --K 4 6 --T 2 3 4 6 8 9 12
  python3 compare_K.py --K 4 6 --dT 0.2 0.8 1.4
  python3 compare_K.py --K 4 6 --x_axis tasks
  python3 compare_K.py --K 4 6 --no_plot
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


def cutoffs_by_K(K_values, base: int, kind: str, scale: str) -> Dict[int, FL.Cutoff]:
    """One cutoff per program number. Under 'genome' scaling the cutoff grows
    in proportion to K, so each program number is read after the same number of
    substitutions per genotype entry; see the module docstring."""
    Ks = sorted(K_values)
    if scale == 'fixed' or kind == 'exposure' or not Ks:
        return {K: FL.Cutoff(kind, base) for K in Ks}
    ref = Ks[0]
    return {K: FL.Cutoff(kind, int(round(base * K / ref))) for K in Ks}


def series(entry, metric: str, dT: float, regime: str, cutoff: FL.Cutoff,
           x_axis: str = 'tasks'):
    """(x, mean, sd) across the T values present for one K, where x is either
    the task number or the task-to-program ratio."""
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
            xs.append(T / spec.K if x_axis == 'ratio' else T)
            ys.append(mu)
            sds.append(sd)
    return np.array(xs), np.array(ys), np.array(sds)


def print_interaction(caches, cutoffs, task_divs):
    """Effect of adding programs, at each task number and in each regime.

    This is the quantity the figure exists to show. A gain that is flat in T
    under sequential selection but grows with T under simultaneous selection is
    the signature of program number limiting one regime and not the other.
    """
    Ks = sorted(caches)
    if len(Ks) < 2:
        return

    print(f'\n{"=" * 96}')
    print('EFFECT OF ADDING PROGRAMS')
    print('Difference in each metric between consecutive program numbers, at '
          'matched task number.')
    print('Flat in T under m=1 and growing in T under m=T is the expected '
          'signature.')
    print('=' * 96)
    head = (f'{"K":>9} {"T":>3} {"dT":>5} '
            f'{"d_diff_m1":>10} {"d_diff_mT":>10} '
            f'{"d_opt_m1":>10} {"d_opt_mT":>10}')
    print(head)
    print('-' * len(head))

    for lo, hi in zip(Ks[:-1], Ks[1:]):
        shared = sorted(set(caches[lo]['spec'].T_values)
                        & set(caches[hi]['spec'].T_values))
        for T in shared:
            for dT in task_divs:
                cells = []
                ok = True
                for metric in ('differentiation', 'optimization'):
                    for sel in ('min', 'T'):
                        vals = []
                        for K in (lo, hi):
                            by_m = caches[K]['data'].get(T, {}).get(dT, {})
                            m = FL.resolve_m(sel, by_m.keys(), T)
                            if m is None:
                                ok = False
                                break
                            _, v = FL.metric_values(by_m[m], metric,
                                                    cutoffs[K])
                            vals.append(FL.mean_sd(v)[0])
                        if not ok:
                            break
                        cells.append(vals[1] - vals[0])
                    if not ok:
                        break
                if not ok or len(cells) < 4:
                    continue
                print(f'{f"{lo}->{hi}":>9} {T:>3} {dT:>5.1f} '
                      f'{cells[0]:>+10.4f} {cells[1]:>+10.4f} '
                      f'{cells[2]:>+10.4f} {cells[3]:>+10.4f}')


def print_table(caches, cutoffs, task_divs):
    print(f'\n{"=" * 100}')
    print('DIFFERENTIATION AND OPTIMIZATION BY PROGRAM AND TASK NUMBER')
    print('Rows at matched T/K are the comparable ones: K=4 at T=8 and K=6 at '
          'T=12 are both two tasks per program.')
    print('cutoff is per program number; see the module docstring.')
    print('=' * 100)
    head = (f'{"K":>3} {"T":>3} {"T/K":>5} {"dT":>5} {"cutoff":>7} '
            f'{"diff_m1":>16} {"diff_mT":>16} {"opt_m1":>16} {"opt_mT":>16}')
    print(head)
    print('-' * len(head))

    for K in sorted(caches):
        spec, data = caches[K]['spec'], caches[K]['data']
        cutoff = cutoffs[K]
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
                      f'{cutoff.value:>7} '
                      f'{cells[0]:>16} {cells[1]:>16} '
                      f'{cells[2]:>16} {cells[3]:>16}')


def make_figure(caches, task_divs, cutoffs, x_axis='ratio', save_path=None):
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
                    xs, ys, sds = series(caches[K], metric, dT, regime,
                                         cutoffs[K], x_axis)
                    if xs.size == 0:
                        continue
                    ax.errorbar(xs, ys, yerr=sds, fmt='o', ls=styles[K],
                                color=colors[dT], lw=0.75, ms=4,
                                markerfacecolor='none', markeredgecolor=colors[dT],
                                capsize=2, capthick=0.8, elinewidth=0.8)

            # T = K, once per program number: to its left that K has at least
            # as many programs as tasks.
            marks = [1.0] if x_axis == 'ratio' else sorted(caches)
            for K in marks:
                ax.axvline(K, color='gray', ls=':', lw=1.0, alpha=0.8, zorder=0)
                if r == 0:
                    lab = '$T=K$' if x_axis == 'ratio' else f'$T=K={K}$'
                    ax.annotate(lab, xy=(K, 1.0),
                                xycoords=('data', 'axes fraction'),
                                xytext=(3, -3), textcoords='offset points',
                                fontsize=8, color='gray', ha='left', va='top')
            ax.axhline(1, color='gray', ls='--', lw=0.8, alpha=0.5)
            ax.set_ylim(0, 1.05)
            ax.set_xlabel('Tasks per program, $T/K$' if x_axis == 'ratio'
                          else 'Number of tasks')
            if x_axis == 'ratio':
                ratios = sorted({T / e['spec'].K for e in caches.values()
                                 for T in e['spec'].T_values})
                ticks = [r for r in (0.5, 1.0, 1.5, 2.0)
                         if min(ratios) - 1e-9 <= r <= max(ratios) + 1e-9]
                if ticks:
                    ax.set_xticks(ticks)
                    ax.set_xticklabels([f'{t:g}' for t in ticks])
            else:
                all_T = sorted({T for e in caches.values()
                                for T in e['spec'].T_values})
                ax.set_xticks(all_T)
                ax.set_xticklabels([str(T) for T in all_T])
            if c == 0:
                ax.set_ylabel(ylabel)
            else:
                ax.tick_params(labelleft=False)

    multi = len({c.value for c in cutoffs.values()}) > 1
    handles = [plt.Line2D([], [], color='0.3', ls=styles[K],
                          label=(f'$K={K}$, {cutoffs[K].value} subs' if multi
                                 else f'$K={K}$'))
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
    p.add_argument('--T', type=int, nargs='+', default=[2, 3, 4, 6, 8, 9, 12],
                   dest='T_values',
                   help='Candidate T values; those without a cache are skipped.')
    p.add_argument('--dT', type=float, nargs='+', default=[0.2, 0.8, 1.4],
                   dest='task_divs')
    p.add_argument('--gamma', type=float, default=1.0)
    p.add_argument('--fitness_r', type=float, default=0.0)
    p.add_argument('--density', type=float, default=0.25)
    p.add_argument('--cutoff', type=int, default=200,
                   help='Cutoff at the smallest program number; scaled per K '
                        'unless --cutoff_scale fixed.')
    p.add_argument('--cutoff_scale', default='genome',
                   choices=['genome', 'fixed'],
                   help="'genome' scales the cutoff with K so every program "
                        "number is read after the same number of substitutions "
                        "per genotype entry; 'fixed' uses one cutoff for all K.")
    p.add_argument('--cutoff_kind', default='substitutions',
                   choices=['substitutions', 'exposure'])
    p.add_argument('--x_axis', default='ratio', choices=['ratio', 'tasks'],
                   help="'ratio' plots against T/K with the boundary at 1 "
                        "(default); 'tasks' plots against T with T=K marked "
                        "per program number.")
    p.add_argument('--no_plot', action='store_true')
    p.add_argument('--no_show', action='store_true')
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    cutoff = FL.Cutoff(args.cutoff_kind, args.cutoff)

    caches = collect(args.K_values, {}, args, cutoff)
    if not caches:
        raise SystemExit('No caches found.')

    cutoffs = cutoffs_by_K(caches.keys(), args.cutoff, args.cutoff_kind,
                           args.cutoff_scale)
    print('\nCutoff per program number: '
          + ', '.join(f'K={K}: {c.label()}' for K, c in sorted(cutoffs.items())))

    print_table(caches, cutoffs, args.task_divs)
    print_interaction(caches, cutoffs, args.task_divs)

    if not args.no_plot:
        os.makedirs(args.save_dir, exist_ok=True)
        path = os.path.join(
            args.save_dir,
            f'{args.filename}_{args.x_axis}_{args.cutoff_scale}'
            f'{args.cutoff}_gamma{args.gamma}_fr{args.fitness_r}'
            f'_density{args.density:.4f}.{args.fmt}')
        fig = make_figure(caches, args.task_divs, cutoffs,
                          x_axis=args.x_axis, save_path=path)
        if args.no_show:
            plt.close(fig)
        else:
            plt.show()