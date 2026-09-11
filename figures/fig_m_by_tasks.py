#!/usr/bin/env python3
"""
fig_m_by_tasks.py
=================
Figure 4: differentiation and optimization against task number, one line per
simultaneity level (m = 1, T/2, T) and one column per task divergence.

A fractional level is plotted only where it gives a whole number of tasks.
Rounding would place a point at a different degree of simultaneity from the one
the line represents, in a figure whose purpose is to hold that degree constant
across task numbers.

Usage:
  python3 fig_m_by_tasks.py --plain_name
  python3 fig_m_by_tasks.py --m 1 1/4 1/2 3/4 T
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
from fractions import Fraction
from typing import List, Optional, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import figlib as FL


METRICS = [('differentiation', 'Degree of differentiation'),
           ('optimization', 'Degree of optimization')]

def m_label(m_spec) -> str:
    if m_spec == 'T':
        return '$m=T$'
    if isinstance(m_spec, Fraction):
        num, den = m_spec.numerator, m_spec.denominator
        return f'$m=T/{den}$' if num == 1 else f'$m={num}T/{den}$'
    return f'$m={m_spec}$'


def m_style(m_spec):
    """Linestyle is the simultaneity channel throughout the paper: dotted at
    the sequential limit, dashed at intermediate m, solid at the simultaneous
    limit. See the visual grammar in figlib."""
    if m_spec == 'T':
        return FL.LS_MT
    if m_spec == 1:
        return FL.LS_M1
    return FL.LS_MHALF


def level_int(m_spec, T: int) -> Optional[int]:
    """The whole number of tasks this level denotes at this T, or None where
    the fraction does not divide the repertoire."""
    if m_spec == 'T':
        return T
    if isinstance(m_spec, Fraction):
        product = m_spec * T
        return int(product) if product.denominator == 1 else None
    return int(m_spec)


def resolve_level(m_spec, available, T: int):
    """Simultaneity level for one task number, or None if the level is not
    defined there.

    A fraction is honoured only when it gives a whole number of tasks. Rounding
    would place a point at a different degree of simultaneity from the one the
    line represents, which is misleading in a figure whose whole purpose is to
    hold that degree constant across task numbers.
    """
    if isinstance(m_spec, Fraction):
        product = m_spec * T
        if product.denominator != 1:
            return None
        target = int(product)
        return target if target in available else None
    return FL.resolve_m(m_spec, available, T)


def make_figure(data, spec: FL.CacheSpec, cutoff: FL.Cutoff,
                m_specs: List, show_K_line: bool = True,
                save_path: Optional[str] = None):
    FL.apply_style()
    n_cols = len(spec.task_divs)
    fig, axes = plt.subplots(2, n_cols, squeeze=False,
                             figsize=(3.4 * n_cols + 0.8, 7.4))
    fig.subplots_adjust(hspace=0.28, wspace=0.18,
                        left=0.10, right=0.97, top=0.90, bottom=0.12)

    # Colour means task divergence here exactly as it does in Figures 2 and 3,
    # even though each column already fixes it: the reader carries one key.
    colors = FL.dt_colors(spec.task_divs)
    t_values = np.array(spec.T_values)

    for r, (metric, ylabel) in enumerate(METRICS):
        for c, dT in enumerate(spec.task_divs):
            ax = axes[r][c]
            ax.set_box_aspect(1)
            ax.text(-0.16, 1.08, FL.panel_label(r * n_cols + c),
                    transform=ax.transAxes, fontsize=14, fontweight='bold',
                    va='top', ha='left')
            if r == 0:
                ax.set_title(fr'$\overline{{\Delta T}}={dT:g}$',
                             fontsize=12, pad=6)

            for i, m_spec in enumerate(m_specs):
                xs, ys, sds = [], [], []
                for T in spec.T_values:
                    by_m = data.get(T, {}).get(dT, {})
                    m = resolve_level(m_spec, by_m.keys(), T)
                    if m is None:
                        continue
                    _, vals = FL.metric_values(by_m[m], metric, cutoff)
                    mu, sd = FL.mean_sd(vals)
                    if np.isfinite(mu):
                        xs.append(T); ys.append(mu); sds.append(sd)
                if not xs:
                    continue
                # Three same-coloured bands overlap in every panel, so the
                # alpha is lower here than in the one-line-per-colour figures.
                FL.band(ax, xs, ys, sds, color=colors[dT],
                        ls=m_style(m_spec), alpha=0.10)

            ax.set_xticks(t_values)
            ax.set_xticklabels([str(int(v)) for v in t_values])
            if r == 1:
                ax.set_xlabel('Number of tasks')

            if show_K_line and t_values.min() <= spec.K <= t_values.max():
                FL.mark_K_line(ax, spec.K, label=f'$K = {spec.K}$',
                               show_label=(r == 0 and c == 0))

            FL.metric_axis(ax, metric, ylabel=False)
            if c == 0:
                ax.set_ylabel(ylabel)
            else:
                ax.tick_params(labelleft=False)

    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color='0.35', ls=m_style(m), lw=1.7,
                      label=m_label(m)) for m in m_specs]
    # One key, in the first panel, next to the one K annotation.
    axes[0][0].legend(handles=handles, fontsize=11, frameon=False,
                      loc='upper right', labelspacing=0.35, handlelength=2.6)

    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
        print(f'Saved: {save_path}')
    return fig


def print_summary(data, spec: FL.CacheSpec, cutoff: FL.Cutoff, m_specs):
    print(f'\n{"=" * 96}')
    print(f'SUMMARY  ({cutoff.label()})  {spec.label()}')
    print('=' * 96)
    head = (f'{"T":>3} {"dT":>5} {"m_spec":>7} {"m":>3} {"n":>5} '
            f'{"differentiation":>18} {"optimization":>18}')
    print(head)
    print('-' * len(head))

    for T in spec.T_values:
        for dT in spec.task_divs:
            by_m = data.get(T, {}).get(dT, {})
            if not by_m:
                continue
            for m_spec in m_specs:
                m = resolve_level(m_spec, by_m.keys(), T)
                if m is None:
                    continue
                _, dv = FL.metric_values(by_m[m], 'differentiation', cutoff)
                _, ov = FL.metric_values(by_m[m], 'optimization', cutoff)
                d_mu, d_sd = FL.mean_sd(dv)
                o_mu, o_sd = FL.mean_sd(ov)
                print(f'{T:>3} {dT:>5.1f} {m_label(m_spec)[3:-1]:>7} '
                      f'{m:>3} {len(dv):>5} '
                      f'{d_mu:>8.4f}+/-{d_sd:<8.4f} '
                      f'{o_mu:>8.4f}+/-{o_sd:<8.4f}')


def parse_m(text: str):
    """'T' is the simultaneous limit, 'a/b' a fraction of the repertoire, and a
    bare integer a fixed number of tasks."""
    if text.upper() == 'T':
        return 'T'
    if '/' in text or '.' in text:
        frac = Fraction(text).limit_denominator(100)
        if not 0 < frac < 1:
            raise ValueError(f'Fraction must lie strictly between 0 and 1: {text}')
        return frac
    return int(text)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--cache_dir', default=_repo_path('simulation_cache'))
    p.add_argument('--save_dir', default=_repo_path('figures_out'))
    p.add_argument('--filename', default='F4')
    p.add_argument('--fmt', default='pdf')
    p.add_argument('--plain_name', action='store_true',
                   help='Write F4.pdf rather than F4_cut200_....pdf, i.e. the '
                        'name the manuscript expects.')
    p.add_argument('--L', type=int, default=100)
    p.add_argument('--K', type=int, default=4)
    p.add_argument('--gamma', type=float, default=1.0)
    p.add_argument('--fitness_r', type=float, default=0.0)
    p.add_argument('--density', type=float, default=0.25)
    p.add_argument('--T', type=int, nargs='+', default=[2, 4, 6, 8],
                   dest='T_values')
    p.add_argument('--dT', type=float, nargs='+', default=[0.2, 0.8, 1.4],
                   dest='task_divs',
                   help='One column per value; a few values keep the panels '
                        'readable.')
    p.add_argument('--m', nargs='+',
                   default=['1', '1/2', 'T'],
                   dest='m_specs',
                   help="Simultaneity levels. '1' and 'T' are the sequential "
                        "and simultaneous limits; 'a/b' is a fraction of the "
                        'task repertoire, plotted only where it gives a whole '
                        'number of tasks; a bare integer is a fixed count.')
    p.add_argument('--cutoff', type=int, default=400)
    p.add_argument('--cutoff_kind', default='substitutions',
                   choices=['substitutions', 'exposure'])
    p.add_argument('--no_K_line', action='store_true')
    p.add_argument('--no_show', action='store_true')
    p.add_argument('--no_summary', action='store_true')
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    m_specs = [parse_m(s) for s in args.m_specs]

    spec = FL.CacheSpec(cache_dir=args.cache_dir, L=args.L, K=args.K,
                        gamma=args.gamma, fitness_r=args.fitness_r,
                        density=args.density, T_values=args.T_values,
                        task_divs=args.task_divs)
    cutoff = FL.Cutoff(args.cutoff_kind, args.cutoff)

    print(f'Loading {spec.label()} ...')
    data = FL.load_grid(
        spec, m_values=lambda T: [level_int(ms, T) for ms in m_specs])

    os.makedirs(args.save_dir, exist_ok=True)
    path = os.path.join(args.save_dir, f'{args.filename}.{args.fmt}') \
        if args.plain_name else os.path.join(
            args.save_dir,
            f'{args.filename}_cut{args.cutoff}_gamma{args.gamma}'
            f'_fr{args.fitness_r}_K{args.K}_density{args.density:.4f}.{args.fmt}')

    fig = make_figure(data, spec, cutoff, m_specs,
                      show_K_line=not args.no_K_line, save_path=path)
    if not args.no_summary:
        print_summary(data, spec, cutoff, m_specs)
    if args.no_show:
        plt.close(fig)
    else:
        plt.show()