#!/usr/bin/env python3
"""
fig_termination.py
==================
Figure S1: where the comparison point falls along the trajectory.

  A, B  the two plotted metrics against substitution number
  C     fraction of single-bit mutants that are beneficial
  D     substitutions to an absorbing state, one line per parameter set

A-C are trajectories at one task number with the two regimes overlaid. Under
simultaneous selection the active set is the whole repertoire, so a genotype
with no beneficial mutation cannot change again: C decays to zero and the
metrics stop moving. Under sequential selection the set is redrawn each epoch
and an absorbing state generally does not exist; C settles at a non-zero level,
meaning the genotype keeps turning over, while A and B level off anyway.
Neither panel establishes that alone: flat metrics by themselves would suggest
convergence to an optimum, which is not what happens.

The two regimes are not on a common scale in C. At m = 1 a mutant need only be
beneficial for the one drawn task; at m = T for all of them jointly. The
contrast to read is the shape, not the height.

Absorbing replicates are held at their last state, which is exact rather than an
extrapolation. Replicates that stopped for any other reason are marked undefined
beyond their last state rather than forward-filled, and a trajectory is cut off
where fewer than MIN_COVERAGE of replicates remain defined.

Bands in A-C are the standard error of the mean, not the spread across
replicates shown in Figures 2-4. The question here is whether the mean
trajectory has stopped moving, which is a statement about the mean; once every
replicate is absorbing the spread across task worlds stops changing and would
stay constant forever, which says nothing about convergence.

D shows how long an absorbing state takes and therefore whether a fixed
comparison point falls before or after it. Everything is read from cached
trajectories and termination records; nothing is re-simulated.

Usage:
  python3 fig_termination.py
  python3 fig_termination.py --T_ref 4 --cutoff 400
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
import glob
import json
import os
import re
from collections import defaultdict
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np

import figlib as FL

SCHEMA_VERSION = 2
ROOT_RE = re.compile(r'L(\d+)_K(\d+)_gamma([\d.]+)_fr(-?[\d.]+)_v(\d+)')


def scan(cache_dir: str, gamma: float, fitness_r: float) -> List[Dict]:
    """One record per condition: parameters plus termination counts and the
    realized substitution counts of the populations that stopped."""
    out = []
    for root in sorted(glob.glob(os.path.join(cache_dir, 'L*_K*_v*'))):
        match = ROOT_RE.search(os.path.basename(root))
        if not match:
            continue
        L, K = int(match.group(1)), int(match.group(2))
        if float(match.group(3)) != gamma or float(match.group(4)) != fitness_r:
            continue

        for dens_dir in sorted(glob.glob(os.path.join(root, 'density*'))):
            rho = float(os.path.basename(dens_dir).replace('density', ''))
            for path in glob.glob(os.path.join(dens_dir, '*_meta.json')):
                try:
                    with open(path) as f:
                        meta = json.load(f)
                except (OSError, ValueError):
                    continue
                if meta.get('schema_version') != SCHEMA_VERSION:
                    continue
                params = meta['params']
                reasons = defaultdict(int)
                stopped = []
                for rep in meta['rep_meta']:
                    reasons[rep['termination_reason']] += 1
                    if rep['termination_reason'] == 'absorbing':
                        stopped.append(rep['n_subs_realized'])
                out.append({
                    'L': L, 'K': K, 'rho': rho,
                    'T': params['T'], 'm': params['m'], 'dT': params['dT'],
                    'n': len(meta['rep_meta']),
                    'n_absorbing': reasons['absorbing'],
                    'subs': stopped,
                })
    return out


def set_label(K: int, rho: float) -> str:
    """Always the numeric density. The K sweep holds rho at 0.25 rather than at
    1/K, so writing 1/K for K=4 would imply a rule the other sets do not
    follow."""
    return fr'$K={K}$, $\rho={rho:g}$' 


TRAJ = {
    'differentiation': 'Degree of differentiation',
    'optimization': 'Degree of optimization',
    'n_ben': 'Fraction of mutations\nthat are beneficial',
}


MIN_COVERAGE = 0.9      # stop a trajectory once this fraction is not defined


def trajectory(reps, kind: str, S: int, LK: int):
    """(values, defined) matrices of one quantity against substitution number.

    An absorbing replicate is held at its last state, which is exact rather
    than an extrapolation: its genotype cannot change again, so that value IS
    its value at every later step. For the beneficial fraction the held value
    is zero, since absorbing means precisely that no beneficial mutation
    exists.

    A replicate that stopped for any other reason -- the substitution budget,
    or the redraw safeguard -- has no defined value beyond its last state, and
    is marked undefined there rather than forward-filled.
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


def panel_traj(ax, data, kind: str, T_ref: int, traj_divs, S: int, LK: int,
               cutoffs=(200,)):
    """One metric against substitution number, sequential over simultaneous."""
    colors = FL.dt_colors(traj_divs)
    x = np.arange(S + 1)

    for dT in traj_divs:
        for m, ls in ((1, FL.LS_M1), (T_ref, FL.LS_MT)):
            reps = FL.get(data, T_ref, dT, m)
            if not reps:
                continue
            M, D = trajectory(reps, kind, S, LK)
            if M.shape[0] == 0:
                continue
            Mn = np.where(D, M, np.nan)
            n = D.sum(axis=0)
            enough = n >= MIN_COVERAGE * D.shape[0]
            with np.errstate(invalid='ignore'):
                mu = np.nanmean(Mn, axis=0)
                se = np.nanstd(Mn, axis=0, ddof=1) / np.sqrt(np.maximum(n, 1))
            mu = np.where(enough, mu, np.nan)
            FL.band(ax, x, mu, se, color=colors[dT], ls=ls, alpha=0.18)

    for c in cutoffs:
        if c <= S:
            FL.mark_cutoff(ax, c)
    ax.set_xlim(0, S)
    ax.set_xlabel('Substitutions')
    ax.set_ylabel(TRAJ[kind])
    if kind in FL.METRIC_YLIM:
        ax.set_ylim(*FL.METRIC_YLIM[kind])
    else:
        ax.set_ylim(bottom=0)


def panel_B(ax, records, sets, cutoff: int, task_divs):
    """Median substitutions to absorption against task divergence, one line per
    parameter set, pooled over task number. The shaded band is the interquartile
    range across replicates."""
    # Task divergence is already the x-axis, so colour is free -- and unused:
    # the parameter sets are told apart by marker alone. Dispersion is drawn as
    # a whisker on each marker rather than as a band, because four bands in one
    # colour could not be attributed to their lines.
    markers = ['o', 's', '^', 'D', 'v']
    dodge = 0.014

    for i, (K, rho) in enumerate(sets):
        xs, meds, los, his = [], [], [], []
        for dT in task_divs:
            vals = [v for r in records
                    if r['K'] == K and abs(r['rho'] - rho) < 1e-6
                    and r['m'] == r['T'] and abs(r['dT'] - dT) < 1e-9
                    for v in r['subs']]
            if not vals:
                continue
            xs.append(dT)
            meds.append(np.median(vals))
            los.append(np.percentile(vals, 25))
            his.append(np.percentile(vals, 75))
        if not xs:
            continue
        off = (i - (len(sets) - 1) / 2.0) * dodge
        xo = [x + off for x in xs]
        ax.vlines(xo, los, his, color='black', lw=0.8, alpha=0.45, zorder=1)
        ax.plot(xo, meds, '-', marker=markers[i % len(markers)], color='black',
                lw=1.0, ms=5, markerfacecolor='none', markeredgecolor='black',
                label=set_label(K, rho), zorder=2)

    ax.set_ylim(top=max(ax.get_ylim()[1], cutoff * 1.15))
    FL.mark_cutoff(ax, cutoff, f'{cutoff} substitutions', axis='y')

    ax.set_xlabel('Mean task divergence')
    ax.set_ylabel('Substitutions to absorbing state')
    ax.set_xticks(task_divs)
    ax.set_xticklabels([f'{v:g}' for v in task_divs])
    ax.legend(fontsize=9, frameon=False, loc='upper left')


def make_figure(records, traj_data, sets, task_divs, traj_divs, T_ref, S, LK,
                cutoff, save_path: Optional[str] = None):
    FL.apply_style()
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 8.0))
    fig.subplots_adjust(wspace=0.30, hspace=0.32, left=0.10, right=0.97,
                        top=0.92, bottom=0.09)
    flat = axes.ravel()

    for i, ax in enumerate(flat):
        ax.text(-0.16, 1.08, FL.panel_label(i), transform=ax.transAxes,
                fontsize=14, fontweight='bold', va='top', ha='left')

    for ax, kind in zip(flat[:3],
                        ('differentiation', 'optimization', 'n_ben')):
        panel_traj(ax, traj_data, kind, T_ref, traj_divs, S, LK,
                   cutoffs=(cutoff,))
    # Both keys sit in the first panel; panel D has its own.
    from matplotlib.lines import Line2D
    cols = FL.dt_colors(traj_divs)
    key_dT = flat[0].legend(
        handles=[Line2D([], [], color=cols[d], lw=1.7,
                        label=fr'$\overline{{\Delta T}} = {d:g}$')
                 for d in traj_divs],
        fontsize=9, frameon=False, loc='upper left')
    flat[0].add_artist(key_dT)
    flat[0].legend(handles=[Line2D([], [], color='0.35', lw=1.7, ls=ls,
                                   label=lab)
                            for ls, lab in ((FL.LS_M1, r'$m = 1$'),
                                            (FL.LS_MT, fr'$m = T = {T_ref}$'))],
                   fontsize=9, frameon=False, loc='upper right',
                   handlelength=2.6)

    panel_B(flat[3], records, sets, cutoff, task_divs)

    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
        print(f'Saved: {save_path}')
    return fig


def print_summary(records, sets, cutoff):
    print(f'\n{"=" * 78}')
    print('TERMINATION SUMMARY')
    print('=' * 78)
    print(f'{"K":>3} {"rho":>6} {"regime":>13} {"n":>7} {"%stopped":>9} '
          f'{"median":>8} {"range":>13}')
    print('-' * 78)

    for K, rho in sets:
        for regime in ('simultaneous', 'sequential'):
            rows = [r for r in records
                    if r['K'] == K and abs(r['rho'] - rho) < 1e-6
                    and ((r['m'] == r['T']) if regime == 'simultaneous'
                         else (r['m'] == 1))]
            if not rows:
                continue
            n = sum(r['n'] for r in rows)
            stopped = sum(r['n_absorbing'] for r in rows)
            subs = [s for r in rows for s in r['subs']]
            med = f'{int(np.median(subs))}' if subs else '--'
            rng = (f'{int(np.min(subs))}-{int(np.max(subs))}' if subs else '--')
            print(f'{K:>3} {rho:>6.3f} {regime:>13} {n:>7} '
                  f'{100.0 * stopped / n:>8.1f}% {med:>8} {rng:>13}')

    print(f'\nMedian substitutions to absorption per condition '
          f'(simultaneous selection), against the {cutoff}-substitution '
          f'comparison:')
    print(f'{"K":>3} {"rho":>6} {"T":>4} {"dT":>5} {"%stopped":>9} {"median":>8}')
    for r in sorted(records, key=lambda r: (r['K'], r['rho'], r['T'], r['dT'])):
        if r['m'] != r['T']:
            continue
        med = int(np.median(r['subs'])) if r['subs'] else 0
        flag = '  <-- after comparison' if med > cutoff else ''
        print(f'{r["K"]:>3} {r["rho"]:>6.3f} {r["T"]:>4} {r["dT"]:>5.1f} '
              f'{100.0 * r["n_absorbing"] / r["n"]:>8.1f}% {med:>8}{flag}')


def parse_set(text: str):
    K = rho = None
    for part in text.split(','):
        key, _, val = part.partition('=')
        if key.strip() == 'K':
            K = int(val)
        elif key.strip() in ('rho', 'density'):
            rho = float(val)
    if K is None:
        raise ValueError(f'Could not parse parameter set: {text!r}')
    return K, (1.0 / K if rho is None else rho)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--cache_dir', default=_repo_path('simulation_cache'))
    p.add_argument('--save_dir', default=_repo_path('figures_out'))
    p.add_argument('--filename', default='FS1_termination')
    p.add_argument('--fmt', default='pdf')
    p.add_argument('--gamma', type=float, default=1.0)
    p.add_argument('--fitness_r', type=float, default=0.0)
    p.add_argument('--K_ref', type=int, default=4,
                   help='Program number shown in the trajectory panels.')
    p.add_argument('--rho_ref', type=float, default=0.25,
                   help='Initial density shown in the trajectory panels.')
    p.add_argument('--T_ref', type=int, default=8,
                   help='Task number for the trajectory panels. The default is '
                        'the largest in the baseline grid, where each phenotype '
                        'is exposed least often under sequential selection and '
                        'the comparison is hardest to defend.')
    p.add_argument('--traj_dT', type=float, nargs='+', default=[0.2, 0.8, 1.4],
                   help='Task divergences drawn in the trajectory panels.')
    p.add_argument('--max_step', type=int, default=800,
                   help='Right-hand limit of the trajectory panels.')
    p.add_argument('--L', type=int, default=100)
    p.add_argument('--sets', nargs='+',
                   default=['K=4,rho=0.25', 'K=4,rho=0.5',
                            'K=6,rho=0.25', 'K=8,rho=0.25'],
                   help='Parameter sets shown in panel B.')
    p.add_argument('--dT', type=float, nargs='+',
                   default=[0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4],
                   dest='task_divs')
    p.add_argument('--cutoff', type=int, default=400)
    p.add_argument('--no_show', action='store_true')
    p.add_argument('--no_summary', action='store_true')
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    sets = [parse_set(s) for s in args.sets]

    print(f'Scanning {args.cache_dir} ...')
    records = scan(args.cache_dir, args.gamma, args.fitness_r)
    if not records:
        raise SystemExit('No conditions found.')
    print(f'  {len(records)} conditions')

    spec = FL.CacheSpec(cache_dir=args.cache_dir, L=args.L, K=args.K_ref,
                        gamma=args.gamma, fitness_r=args.fitness_r,
                        density=args.rho_ref, T_values=[args.T_ref],
                        task_divs=args.traj_dT)
    print(f'Loading trajectories for T={args.T_ref} ...')
    traj_data = FL.load_grid(spec, m_values=lambda T: [1, T], verbose=False)

    os.makedirs(args.save_dir, exist_ok=True)
    path = os.path.join(args.save_dir, f'{args.filename}.{args.fmt}')
    fig = make_figure(records, traj_data, sets, args.task_divs, args.traj_dT,
                      args.T_ref, args.max_step, args.L * args.K_ref,
                      args.cutoff, save_path=path)

    if not args.no_summary:
        print_summary(records, sets, args.cutoff)
    if args.no_show:
        plt.close(fig)
    else:
        plt.show()