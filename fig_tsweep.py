#!/usr/bin/env python3
"""
fig_tsweep.py
=============
Differentiation and optimization against task number, at two points along the
trajectory, with the simultaneity m held fixed across the whole figure.

F2 and F3 differ only in that setting: F2 is sequential selection (m=1) and F3
is simultaneous selection (m=T). They are generated from this one module so
that the two panels a reader is asked to compare cannot diverge in how a cutoff
is resolved or how differentiation is normalized.

  Row 1  degree of differentiation, mean pairwise phenotype distance over each
         replicate's own realized task divergence
  Row 2  degree of optimization, 1 - ||d||_2 / sqrt(T)
  Cols   two cutoffs, early and late
  Lines  one per task divergence, viridis_r

Points are means across replicates and bars are +/- 1 SD. Each replicate is an
independent draw of an initial genome and a task ensemble, so the SD describes
variation in the plotted quantity across task worlds rather than path noise
within one.

The vertical marker at T = K matters for the argument: to its left there are at
least as many programs as tasks, so a decline there cannot be explained by a
shortage of programs.

Usage:
  python3 fig_tsweep.py --which F2
  python3 fig_tsweep.py --which F3
  python3 fig_tsweep.py --which F2 --cutoff_kind exposure --cutoffs 50 100
  python3 fig_tsweep.py --which F3 --gamma 4.0 --filename FS_gamma4_mT
"""

import argparse
import os
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

import figlib as FL


# ============================================================
# CONFIG
# ============================================================

@dataclass
class FigConfig:
    m_selector: object = 'min'          # 'min' = sequential, 'T' = simultaneous
    figsize: Tuple[float, float] = (9.0, 8.0)
    phase_labels: List[str] = field(
        default_factory=lambda: ['Early Phase', 'Late Phase'])
    show_K_line: bool = True
    K_line_label: bool = True
    annotate_cutoff: bool = True


PRESETS = {
    'F2': dict(m_selector='min', filename='F2',
               phase_labels=['Early Phase', 'Late Phase']),
    'F3': dict(m_selector='T', filename='F3',
               phase_labels=['Early Phase', 'Late Phase']),
}


# ============================================================
# FIGURE
# ============================================================

def make_figure(data, spec: FL.CacheSpec, cutoffs: List[FL.Cutoff],
                fig_cfg: FigConfig, save_path: Optional[str] = None):
    FL.apply_style()
    fig, axes = plt.subplots(2, len(cutoffs), figsize=fig_cfg.figsize,
                             squeeze=False)
    fig.subplots_adjust(hspace=0.30, wspace=0.22,
                        left=0.12, right=0.94, top=0.90, bottom=0.14)

    colors = FL.dt_colors(spec.task_divs)
    t_values = np.array(spec.T_values)
    rows = ['differentiation', 'optimization']

    for col, cutoff in enumerate(cutoffs):
        for row, metric in enumerate(rows):
            ax = axes[row][col]
            ax.set_box_aspect(1)
            ax.text(-0.15, 1.08, FL.panel_label(row * len(cutoffs) + col),
                    transform=ax.transAxes, fontsize=14, fontweight='bold',
                    va='top', ha='left')

            if row == 0:
                title = (fig_cfg.phase_labels[col]
                         if col < len(fig_cfg.phase_labels) else cutoff.label())
                if fig_cfg.annotate_cutoff:
                    title = f'{title}\n({cutoff.label()})'
                ax.set_title(title, fontsize=12, pad=6)

            for dT in spec.task_divs:
                xs, ys, sds = [], [], []
                for T in spec.T_values:
                    by_m = data.get(T, {}).get(dT, {})
                    m = FL.resolve_m(fig_cfg.m_selector, by_m.keys(), T)
                    if m is None:
                        continue
                    _, vals = FL.metric_values(by_m[m], metric, cutoff)
                    mu, sd = FL.mean_sd(vals)
                    if np.isfinite(mu):
                        xs.append(T); ys.append(mu); sds.append(sd)

                if xs:
                    ax.errorbar(xs, ys, yerr=sds, fmt='-o', color=colors[dT],
                                lw=0.75, ms=5, markerfacecolor='none',
                                markeredgecolor=colors[dT], capsize=2,
                                capthick=1.0, elinewidth=1.0)

            ax.set_xticks(t_values)
            ax.set_xticklabels([str(int(v)) for v in t_values])
            ax.set_xlabel('Number of tasks')

            if fig_cfg.show_K_line and t_values.min() <= spec.K <= t_values.max():
                ax.axvline(spec.K, color='gray', ls=':', lw=1.0, alpha=0.8,
                           zorder=0)
                if fig_cfg.K_line_label and row == 0:
                    ax.annotate(f'$T=K={spec.K}$', xy=(spec.K, 1.0),
                                xycoords=('data', 'axes fraction'),
                                xytext=(3, -3), textcoords='offset points',
                                fontsize=8, color='gray', ha='left', va='top')

            ax.axhline(1, color='gray', ls='--', lw=0.8, alpha=0.5)
            ax.set_ylim(0, 1.05)
            if col == 0:
                ax.set_ylabel(FL.metric_label(metric))
            else:
                ax.tick_params(labelleft=False)

    FL.add_dt_colorbar(fig, spec.task_divs)
    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
        print(f'Saved: {save_path}')
    return fig


# ============================================================
# SUMMARY
# ============================================================

def print_summary(data, spec: FL.CacheSpec, cutoffs, fig_cfg: FigConfig):
    print(f'\n{"=" * 118}')
    print(f'SUMMARY  m_selector={fig_cfg.m_selector!r}  {spec.label()}')
    print('dT_real = realized task divergence (mean +/- SD over replicates); '
          'ep = selective epochs per task at the cutoff')
    print('clamp = replicates read at their last state; exact for absorbing '
          'terminations, reportable only when cap > 0')
    print('=' * 118)
    head = (f'{"T":>3} {"dT":>5} {"m":>3} {"cutoff":>26} {"n":>4} {"clamp":>5} '
            f'{"cap":>4} {"dT_real":>15} {"ep_mn":>6} {"ep_min":>6} '
            f'{"diff":>16} {"opt":>16}')
    print(head)
    print('-' * len(head))

    for T in spec.T_values:
        for dT in spec.task_divs:
            by_m = data.get(T, {}).get(dT, {})
            m = FL.resolve_m(fig_cfg.m_selector, by_m.keys(), T)
            if m is None:
                continue
            reps = by_m[m]
            td_mu, td_sd = FL.realized_dT_summary(reps)

            for cutoff in cutoffs:
                term = FL.termination_summary(reps, cutoff)
                expo = FL.exposure_summary(reps, cutoff)
                _, dv = FL.metric_values(reps, 'differentiation', cutoff)
                _, ov = FL.metric_values(reps, 'optimization', cutoff)
                d_mu, d_sd = FL.mean_sd(dv)
                o_mu, o_sd = FL.mean_sd(ov)
                print(f'{T:>3} {dT:>5.1f} {m:>3} {cutoff.label():>26} '
                      f'{term["n_reps"]:>4} {term["n_clamped"]:>5} '
                      f'{term["n_redraw_cap"]:>4} '
                      f'{td_mu:>7.4f}+/-{td_sd:<6.4f} '
                      f'{expo["ep_mean"]:>6.1f} {expo["ep_min"]:>6.1f} '
                      f'{d_mu:>7.4f}+/-{d_sd:<7.4f} '
                      f'{o_mu:>7.4f}+/-{o_sd:<7.4f}')


# ============================================================
# CLI
# ============================================================

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--which', default='F2', choices=list(PRESETS),
                   help='F2 = sequential (m=1), F3 = simultaneous (m=T).')
    p.add_argument('--cache_dir', default='simulation_cache')
    p.add_argument('--save_dir', default='figures_out')
    p.add_argument('--filename', default=None)
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
    p.add_argument('--cutoffs', type=int, nargs='+', default=[50, 200])
    p.add_argument('--cutoff_kind', default='substitutions',
                   choices=['substitutions', 'exposure'],
                   help='substitutions: same evolutionary change (primary). '
                        'exposure: same selective epochs per task (control).')
    p.add_argument('--exposure_mode', default='realized',
                   choices=['realized', 'expected'])
    p.add_argument('--no_show', action='store_true')
    p.add_argument('--no_summary', action='store_true')
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    preset = PRESETS[args.which]

    spec = FL.CacheSpec(cache_dir=args.cache_dir, L=args.L, K=args.K,
                        gamma=args.gamma, fitness_r=args.fitness_r,
                        density=args.density, T_values=args.T_values,
                        task_divs=args.task_divs)
    fig_cfg = FigConfig(m_selector=preset['m_selector'],
                        phase_labels=preset['phase_labels'])
    cutoffs = [FL.Cutoff(args.cutoff_kind, c, args.exposure_mode)
               for c in args.cutoffs]

    print(f'Loading {spec.label()} ...')
    want = None if args.which == 'F2' else None      # load all m, pick later
    data = FL.load_grid(spec, m_values=want)

    os.makedirs(args.save_dir, exist_ok=True)
    stem = args.filename or preset['filename']
    tag = ('' if args.cutoff_kind == 'substitutions'
           else f'_exposure{args.exposure_mode}')
    path = os.path.join(
        args.save_dir,
        f'{stem}{tag}_gamma{args.gamma}_fr{args.fitness_r}'
        f'_K{args.K}_density{args.density:.4f}.{args.fmt}')

    fig = make_figure(data, spec, cutoffs, fig_cfg, save_path=path)
    if not args.no_summary:
        print_summary(data, spec, cutoffs, fig_cfg)
    if args.no_show:
        plt.close(fig)
    else:
        plt.show()
