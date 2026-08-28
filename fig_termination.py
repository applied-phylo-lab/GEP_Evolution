#!/usr/bin/env python3
"""
fig_termination.py
==================
When do simulations stop, and does the comparison point fall after that?

  Panel A  fraction of populations reaching an absorbing state, against task
           number, separately for sequential and simultaneous selection
  Panel B  substitutions required to reach that state under simultaneous
           selection, against task divergence, one line per parameter set,
           with the comparison point marked

Panel A shows that the two selection regimes differ in whether adaptation ends
at all. Under simultaneous selection the active set is the whole task
repertoire, so a genotype with no beneficial mutation cannot change again.
Under sequential selection the active set is redrawn each epoch, and a mutation
that is not beneficial for the currently drawn task may be beneficial for
another, so an endpoint in that sense does not generally exist.

Panel B shows how long reaching that state takes, and therefore whether a
comparison at a fixed number of substitutions falls before or after it. The
requirement depends jointly on the genotype matrix and on how separated the
task optima are: weakly separated optima are satisfied quickly under every
parameter set, whereas strongly separated optima require substantially more
substitutions when the genotype matrix is larger or denser. A comparison point
adequate for one parameter set therefore need not be adequate for another, and
the shortfall falls on the strongly separated conditions.

Both panels are read from the termination records stored with each replicate;
nothing is re-simulated.

Usage:
  python3 fig_termination.py
  python3 fig_termination.py --cutoff 200
  python3 fig_termination.py --sets "K=4,rho=0.25" "K=6,rho=0.25"
"""

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
    return f'$K={K}$, ' + (r'$\rho=1/K$' if abs(rho - 1.0 / K) < 1e-6
                           else fr'$\rho={rho:g}$')


def panel_A(ax, records, K_ref: int, rho_ref: float, task_divs):
    """Fraction stopping, against task number, for the reference parameter set."""
    colors = FL.dt_colors(task_divs)
    styles = {'sequential': ':', 'simultaneous': '-'}

    for regime, marker in (('simultaneous', 'o'), ('sequential', 's')):
        for dT in task_divs:
            rows = [r for r in records
                    if r['K'] == K_ref and abs(r['rho'] - rho_ref) < 1e-6
                    and abs(r['dT'] - dT) < 1e-9
                    and ((r['m'] == r['T']) if regime == 'simultaneous'
                         else (r['m'] == 1))]
            if not rows:
                continue
            rows.sort(key=lambda r: r['T'])
            xs = [r['T'] for r in rows]
            ys = [100.0 * r['n_absorbing'] / r['n'] for r in rows]
            ax.plot(xs, ys, styles[regime], marker=marker, color=colors[dT],
                    lw=0.9, ms=4, markerfacecolor='none',
                    markeredgecolor=colors[dT])

    ax.set_xlabel('Number of tasks')
    ax.set_ylabel('Populations reaching an\nabsorbing state (%)')
    ax.set_ylim(-4, 104)
    Ts = sorted({r['T'] for r in records
                 if r['K'] == K_ref and abs(r['rho'] - rho_ref) < 1e-6})
    ax.set_xticks(Ts)
    ax.set_xticklabels([str(T) for T in Ts])

    handles = [
        plt.Line2D([], [], color='0.3', ls='-', marker='o',
                   markerfacecolor='none', label='Simultaneous ($m=T$)'),
        plt.Line2D([], [], color='0.3', ls=':', marker='s',
                   markerfacecolor='none', label='Sequential ($m=1$)'),
    ]
    ax.legend(handles=handles, fontsize=8, frameon=False, loc='center right')


def panel_B(ax, records, sets, cutoff: int, task_divs):
    """Median substitutions to absorption against task divergence, one line per
    parameter set, pooled over task number. The shaded band is the interquartile
    range across replicates."""
    markers = ['o', 's', '^', 'D', 'v']
    cmap = plt.get_cmap('viridis')

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
        color = cmap(0.08 + 0.78 * i / max(len(sets) - 1, 1))
        ax.fill_between(xs, los, his, color=color, alpha=0.15, linewidth=0)
        ax.plot(xs, meds, '-', marker=markers[i % len(markers)], color=color,
                lw=1.1, ms=4.5, markerfacecolor='none', markeredgecolor=color,
                label=set_label(K, rho))

    ax.axhline(cutoff, color='firebrick', ls='--', lw=1.0, zorder=0)
    ax.annotate(f'{cutoff} substitutions', xy=(0.02, cutoff),
                xycoords=('axes fraction', 'data'),
                xytext=(0, 4), textcoords='offset points',
                fontsize=8, color='firebrick', ha='left', va='bottom')

    ax.set_xlabel('Mean task divergence')
    ax.set_ylabel('Substitutions to absorbing state')
    ax.set_xticks(task_divs)
    ax.set_xticklabels([f'{v:g}' for v in task_divs])
    ax.legend(fontsize=8, frameon=False, loc='upper left')


def make_figure(records, K_ref, rho_ref, sets, task_divs, cutoff,
                save_path: Optional[str] = None):
    FL.apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2))
    fig.subplots_adjust(wspace=0.30, left=0.09, right=0.97,
                        top=0.90, bottom=0.16)

    for i, ax in enumerate(axes):
        ax.text(-0.14, 1.08, FL.panel_label(i), transform=ax.transAxes,
                fontsize=14, fontweight='bold', va='top', ha='left')

    panel_A(axes[0], records, K_ref, rho_ref, task_divs)
    panel_B(axes[1], records, sets, cutoff, task_divs)

    FL.add_dt_colorbar(fig, task_divs, rect=[0.09, 0.005, 0.34, 0.020])

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
    p.add_argument('--cache_dir', default='simulation_cache')
    p.add_argument('--save_dir', default='figures_out')
    p.add_argument('--filename', default='FS1_termination')
    p.add_argument('--fmt', default='pdf')
    p.add_argument('--gamma', type=float, default=1.0)
    p.add_argument('--fitness_r', type=float, default=0.0)
    p.add_argument('--K_ref', type=int, default=4,
                   help='Program number shown in panel A.')
    p.add_argument('--rho_ref', type=float, default=0.25,
                   help='Initial density shown in panel A.')
    p.add_argument('--sets', nargs='+',
                   default=['K=4,rho=0.25', 'K=4,rho=0.5',
                            'K=6,rho=0.25', 'K=8,rho=0.25'],
                   help='Parameter sets shown in panel B.')
    p.add_argument('--dT', type=float, nargs='+',
                   default=[0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4],
                   dest='task_divs')
    p.add_argument('--cutoff', type=int, default=200)
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

    os.makedirs(args.save_dir, exist_ok=True)
    path = os.path.join(args.save_dir, f'{args.filename}.{args.fmt}')
    fig = make_figure(records, args.K_ref, args.rho_ref, sets,
                      args.task_divs, args.cutoff, save_path=path)

    if not args.no_summary:
        print_summary(records, sets, args.cutoff)
    if args.no_show:
        plt.close(fig)
    else:
        plt.show()