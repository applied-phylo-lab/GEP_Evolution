#!/usr/bin/env python3
"""
fig_termination.py
==================
Figure S1: phenotype and beneficial-mutation trajectories at m = 1, T/2, T.

A-C show differentiation, optimization and the beneficial-mutant fraction.
Colour identifies task separation; dotted, dashed and solid lines identify
sequential, intermediate and simultaneous selection. No simulations are run.

Only absorbing endpoints are carried forward. Budget and safeguard stops are
undefined beyond their recorded data, and curves end when fewer than 90% of
replicates remain defined. The beneficial fraction is recorded for the
successful active set preceding an accepted substitution.
"""
import argparse
import os
import sys
from pathlib import Path

_REPO_ROOT = str(Path(__file__).resolve().parents[1])
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

def _repo_path(*parts):
    return os.path.join(_REPO_ROOT, *parts)

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import figlib as FL

TRAJ = {
    'differentiation': 'Degree of differentiation',
    'optimization': 'Degree of optimization',
    'n_ben': 'Fraction of mutations\nthat are beneficial',
}


MIN_COVERAGE = 0.9      # stop a trajectory once this fraction is not defined


def trajectory(reps, kind: str, S: int, LK: int):
    """(values, defined) matrices of one quantity against substitution number.

        An absorbing replicate is held at its last state, which is exact: its
        genotype cannot change again. For the beneficial fraction the held value is
        zero. A replicate that stopped for any other reason is marked undefined
        beyond its last state rather than forward-filled.
    """
    rows, defined = [], []
    for rep in reps:
        if kind == 'differentiation':
            td = float(rep['task_dT_realized'])
            if not (np.isfinite(td) and td > 0):
                continue
            v = np.asarray(rep['pheno_dist'], dtype=float) / td
        elif kind == 'optimization':
            d = np.asarray(rep['d'], dtype=float)
            v = 1.0 - np.linalg.norm(d, axis=1) / np.sqrt(d.shape[1])
        else:
            v = np.asarray(rep['n_ben'], dtype=float) / float(LK)
            if rep.get('termination_reason') == 'absorbing':
                v = np.append(v, 0.0)
        n = v.shape[0]
        if n == 0:
            continue
        steps = np.arange(S + 1)
        rows.append(v[np.minimum(steps, n - 1)])
        defined.append(np.ones(S + 1, dtype=bool)
                       if rep.get('termination_reason') == 'absorbing'
                       else steps < n)
    if not rows:
        return np.empty((0, S + 1)), np.empty((0, S + 1), dtype=bool)
    return np.vstack(rows), np.vstack(defined)


def regimes(T):
    """Preserve the endpoints; add exact half-simultaneity for even T > 2."""
    if T < 2 or T % 2:
        raise ValueError('T_ref must be even and at least 2 for m=T/2.')
    result = [(1, FL.LS_M1, r'$m=1$')]
    if T // 2 != 1:
        result.append((T // 2, FL.LS_MHALF, fr'$m=T/2={T // 2}$'))
    result.append((T, FL.LS_MT, fr'$m=T={T}$'))
    return result


def summarize_trajectory(reps, kind, S, LK):
    values, defined = trajectory(reps, kind, S, LK)
    if values.shape[0] != len(reps):
        raise ValueError('A replicate has missing or invalid trajectory data.')
    n = defined.sum(axis=0)
    total = np.where(defined, values, 0.0).sum(axis=0)
    mean = np.divide(total, n, out=np.full(S + 1, np.nan), where=n > 0)
    squared = np.where(defined, (values - mean) ** 2, 0.0).sum(axis=0)
    var = np.divide(squared, n - 1, out=np.full(S + 1, np.nan), where=n > 1)
    se = np.sqrt(np.divide(var, n, out=np.full(S + 1, np.nan), where=n > 0))
    enough = n >= MIN_COVERAGE * len(reps)
    return np.where(enough, mean, np.nan), np.where(enough, se, np.nan), n


def validate_grid(data, T, divergences):
    for dt in divergences:
        for m, _, _ in regimes(T):
            reps = FL.get(data, T, dt, m)
            if not reps or len(reps) != 200:
                raise ValueError(f'Expected 200 replicates at T={T}, dT={dt}, m={m}.')


def make_figure(data, divergences, T, S, LK, cutoff, early_cutoff=50,
                save_path=None):
    validate_grid(data, T, divergences)
    FL.apply_style()
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0), sharex=True)
    fig.subplots_adjust(left=0.065, right=0.985, bottom=0.18, top=0.88, wspace=0.36)
    colors = FL.dt_colors(divergences)
    x = np.arange(S + 1)
    for i, (ax, kind) in enumerate(zip(axes, TRAJ)):
        for dt in divergences:
            for m, ls, _ in regimes(T):
                mu, se, _ = summarize_trajectory(FL.get(data, T, dt, m), kind, S, LK)
                FL.band(ax, x, mu, se, color=colors[dt], ls=ls, alpha=0.12)
        for c in sorted(set((early_cutoff, cutoff))):
            if 0 <= c <= S:
                ax.axvline(c, color=FL.CUTOFF_COLOR, lw=1.1, alpha=0.8)
                if i == 0:
                    # label inboard of each line, so neither runs into the
                    # legend in the upper right
                    side = -1 if c >= 0.5 * S else 1
                    ax.text(c + side * 0.012 * S, 0.985, f'{c} substitutions',
                            transform=ax.get_xaxis_transform(),
                            rotation=90, ha='left' if side > 0 else 'right',
                            va='top', fontsize=8, color=FL.CUTOFF_COLOR)
        ax.set_xlim(0, S)
        ax.set_xlabel('Substitutions')
        ax.set_xticks(np.arange(0, S + 1, 200))
        ax.set_ylabel(TRAJ[kind])
        ax.set_ylim(*FL.METRIC_YLIM[kind]) if kind in FL.METRIC_YLIM else ax.set_ylim(bottom=0)
        ax.text(-0.13, 1.06, FL.panel_label(i), transform=ax.transAxes,
                fontsize=14, fontweight='bold', va='bottom')
    handles = [Line2D([], [], color=colors[d], lw=1.7,
                     label=fr'$\overline{{\Delta T}}={d:g}$')
               for d in divergences]
    handles.extend(Line2D([], [], color='0.25', lw=1.7, ls=ls, label=lab)
                   for _, ls, lab in regimes(T))
    axes[0].legend(handles=handles, loc='upper right', ncol=1,
                   fontsize=8.5, handlelength=2.2, labelspacing=0.32,
                   borderpad=0.3, frameon=True, facecolor='white',
                   edgecolor='none', framealpha=1.0)
    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
        print(f'Saved: {save_path}')
    return fig


def print_summary(data, divergences, T, S, LK, cutoff):
    print('Recorded trajectory coverage (no forward-filling of budget/safeguard stops)')
    for dt in divergences:
        for m, _, _ in regimes(T):
            for kind in TRAJ:
                mu, _, n = summarize_trajectory(FL.get(data, T, dt, m), kind, S, LK)
                defined = np.flatnonzero(np.isfinite(mu))
                last = int(defined[-1]) if len(defined) else -1
                print(f'dT={dt:g} m={m} {kind}: last plotted={last}, '
                      f'n at cutoff={n[min(cutoff, S)]}')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache_dir', default=_repo_path('simulation_cache'))
    p.add_argument('--save_dir', default=_repo_path('figures_out'))
    p.add_argument('--filename', default='FS1_termination')
    p.add_argument('--fmt', default='pdf')
    p.add_argument('--gamma', type=float, default=1.0)
    p.add_argument('--fitness_r', type=float, default=0.0)
    p.add_argument('--K_ref', type=int, default=4)
    p.add_argument('--rho_ref', type=float, default=0.25)
    p.add_argument('--T_ref', type=int, default=8)
    p.add_argument('--traj_dT', type=float, nargs='+', default=[0.2, 0.8, 1.4])
    p.add_argument('--max_step', type=int, default=800)
    p.add_argument('--L', type=int, default=100)
    p.add_argument('--cutoff', type=int, default=400)
    p.add_argument('--early_cutoff', type=int, default=50)
    p.add_argument('--no_show', action='store_true')
    p.add_argument('--no_summary', action='store_true')
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    spec = FL.CacheSpec(cache_dir=args.cache_dir, L=args.L, K=args.K_ref,
                        gamma=args.gamma, fitness_r=args.fitness_r,
                        density=args.rho_ref, T_values=[args.T_ref],
                        task_divs=args.traj_dT)
    data = FL.load_grid(spec, m_values=[m for m, _, _ in regimes(args.T_ref)], verbose=False)
    os.makedirs(args.save_dir, exist_ok=True)
    path = os.path.join(args.save_dir, f'{args.filename}.{args.fmt}')
    fig = make_figure(data, args.traj_dT, args.T_ref, args.max_step,
                      args.L * args.K_ref, args.cutoff, args.early_cutoff, path)
    if not args.no_summary:
        print_summary(data, args.traj_dT, args.T_ref, args.max_step,
                      args.L * args.K_ref, args.cutoff)
    if args.no_show:
        plt.close(fig)
    else:
        plt.show()
