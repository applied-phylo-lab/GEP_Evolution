#!/usr/bin/env python3
"""
fig_regimes.py
==============
Figures S2, S3 and S5: sequential against simultaneous selection at a single
cutoff. Rows are the two metrics, columns the two regimes, lines one per task
divergence.

This is the layout for every supplementary robustness analysis, so each is one
command producing one file. The main-text figures instead spend their columns on
the early and late phases, because the change between them is part of the
result; a robustness check needs only the endpoint.

Usage:
  python3 fig_regimes.py --fitness_r -2.0     # S2
  python3 fig_regimes.py --gamma 4.0          # S3
  python3 fig_regimes.py --density 0.5        # S5
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
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

import figlib as FL


REGIMES = [('min', 'Sequential ($m=1$)'),
           ('T', 'Simultaneous ($m=T$)')]

METRICS = [('differentiation', 'Degree of differentiation'),
           ('optimization', 'Degree of optimization')]


def make_figure(data, spec: FL.CacheSpec, cutoff: FL.Cutoff,
                show_K_line: bool = True,
                save_path: Optional[str] = None):
    FL.apply_style()
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 8.0), squeeze=False)
    fig.subplots_adjust(hspace=0.30, wspace=0.22,
                        left=0.12, right=0.94, top=0.90, bottom=0.14)

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
                if xs:
                    # One m per panel, so every line is solid; see the visual
                    # grammar in figlib.
                    FL.band(ax, xs, ys, sds, color=colors[dT], ls=FL.LS_MT)

            ax.set_xticks(t_values)
            ax.set_xticklabels([str(int(v)) for v in t_values])
            ax.set_xlabel('Number of tasks')

            if show_K_line and t_values.min() <= spec.K <= t_values.max():
                FL.mark_K_line(ax, spec.K, label=f'$K = {spec.K}$',
                               show_label=(r == 0 and c == 0))

            FL.metric_axis(ax, metric, ylabel=False)
            if c == 0:
                ax.set_ylabel(ylabel)
            else:
                ax.tick_params(labelleft=False)

    FL.add_dt_colorbar(fig, spec.task_divs)
    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
        print(f'Saved: {save_path}')
    return fig


def print_summary(data, spec: FL.CacheSpec, cutoff: FL.Cutoff):
    print(f'\n{"=" * 112}')
    print(f'SUMMARY  sequential vs simultaneous  ({cutoff.label()})  '
          f'{spec.label()}')
    print('dT_real = realized task divergence (mean +/- SD across replicates)')
    print('cap = replicates ended by the redraw safeguard; should be zero')
    print('=' * 112)
    head = (f'{"T":>3} {"dT":>5} {"n":>5} {"cap":>4} {"dT_real":>15} '
            f'{"diff_m1":>16} {"diff_mT":>16} '
            f'{"opt_m1":>16} {"opt_mT":>16}')
    print(head)
    print('-' * len(head))

    for T in spec.T_values:
        for dT in spec.task_divs:
            by_m = data.get(T, {}).get(dT, {})
            if not by_m:
                continue
            any_reps = next(iter(by_m.values()))
            td_mu, td_sd = FL.realized_dT_summary(any_reps)
            cap = sum(FL.termination_summary(reps)['n_redraw_cap']
                      for reps in by_m.values())

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

            print(f'{T:>3} {dT:>5.1f} {len(any_reps):>5} {cap:>4} '
                  f'{td_mu:>7.4f}+/-{td_sd:<6.4f} '
                  f'{cells[0]:>16} {cells[1]:>16} '
                  f'{cells[2]:>16} {cells[3]:>16}')


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--cache_dir', default=_repo_path('simulation_cache'))
    p.add_argument('--save_dir', default=_repo_path('figures_out'))
    p.add_argument('--filename', default='FS')
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
    p.add_argument('--cutoff', type=int, default=400,
                   help='Substitutions, or selective epochs per task when '
                        '--cutoff_kind is exposure.')
    p.add_argument('--cutoff_kind', default='substitutions',
                   choices=['substitutions', 'exposure'])
    p.add_argument('--exposure_mode', default='realized',
                   choices=['realized', 'expected'])
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
    cutoff = FL.Cutoff(args.cutoff_kind, args.cutoff, args.exposure_mode)

    print(f'Loading {spec.label()} ...')
    # Only the two limits are plotted; loading the intermediate levels was
    # enough to exhaust memory at T = 8.
    data = FL.load_grid(spec, m_values=lambda T: [1, T])

    os.makedirs(args.save_dir, exist_ok=True)
    tag = ('' if args.cutoff_kind == 'substitutions'
           else f'_exposure{args.exposure_mode}')
    path = os.path.join(
        args.save_dir,
        f'{args.filename}{tag}_cut{args.cutoff}'
        f'_gamma{args.gamma}_fr{args.fitness_r}'
        f'_K{args.K}_density{args.density:.4f}.{args.fmt}')

    fig = make_figure(data, spec, cutoff,
                      show_K_line=not args.no_K_line, save_path=path)
    if not args.no_summary:
        print_summary(data, spec, cutoff)
    if args.no_show:
        plt.close(fig)
    else:
        plt.show()