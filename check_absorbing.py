#!/usr/bin/env python3
"""
check_absorbing.py
==================
How often did simulations end because no further mutation was beneficial?

Reads only the `_meta.json` sidecars, so it is fast and touches no trajectory
data. It walks every parameter root present in the cache, so it covers program
numbers and task grids that were run outside the named batch specifications.

The distinction it reports is the one the late-phase comparison rests on. Under
simultaneous selection (m = T) the active task set is the whole repertoire, so
a genotype with no beneficial mutation cannot change again and the trajectory
has an endpoint. Under sequential and partially simultaneous selection (m < T)
the active set is redrawn each epoch, so a mutation that is not beneficial for
the currently drawn tasks may be beneficial for another draw; an endpoint in
that sense generally does not exist.

Termination reasons:
  absorbing    no mutation beneficial under any admissible active task set
  n_subs       the substitution budget was reached
  redraw_cap   numerical safeguard; should not appear

Read the `%abs` column for m = T. Values near 100 support describing the late
phase as falling beyond the point at which adaptation has run its course.
Substantially lower values mean the budget was still binding for some
populations, and the claim should be qualified accordingly.

Usage:
  python3 check_absorbing.py
  python3 check_absorbing.py --regime sequential
  python3 check_absorbing.py --dT 1.4
  python3 check_absorbing.py --cache_dir simulation_cache
"""

import argparse
import glob
import json
import os
from collections import defaultdict
from typing import Dict, Optional

SCHEMA_VERSION = 2


def load_meta(path: str) -> Optional[dict]:
    try:
        with open(path) as f:
            meta = json.load(f)
    except (OSError, ValueError):
        return None
    if meta.get('schema_version') != SCHEMA_VERSION:
        return None
    if 'params' not in meta or 'rep_meta' not in meta:
        return None
    return meta


def keep(params: dict, regime: str) -> bool:
    if regime == 'simultaneous':
        return params['m'] == params['T']
    if regime == 'sequential':
        return params['m'] == 1
    return True


def report(cache_dir: str, regime: str, only_dT: Optional[float],
           show_subs: bool):
    roots = sorted(glob.glob(os.path.join(cache_dir, 'L*_K*_v2')))
    if not roots:
        print(f'No parameter roots found under {cache_dir}/')
        return

    grand: Dict[str, int] = defaultdict(int)

    for root in roots:
        for dens in sorted(glob.glob(os.path.join(root, 'density*'))):
            rows: Dict[tuple, Dict[str, int]] = defaultdict(
                lambda: defaultdict(int))
            subs: Dict[tuple, list] = defaultdict(list)

            for path in glob.glob(os.path.join(dens, '*_meta.json')):
                meta = load_meta(path)
                if meta is None:
                    continue
                params = meta['params']
                if not keep(params, regime):
                    continue
                if only_dT is not None and abs(params['dT'] - only_dT) > 1e-9:
                    continue
                key = (params['T'], params['m'], params['dT'])
                for rep in meta['rep_meta']:
                    rows[key][rep['termination_reason']] += 1
                    grand[rep['termination_reason']] += 1
                    subs[key].append(rep['n_subs_realized'])

            if not rows:
                continue

            print(f'\n{os.path.basename(root)}  {os.path.basename(dens)}  '
                  f'[{regime}]')
            head = (f'{"T":>3} {"m":>3} {"dT":>5} {"n":>5} '
                    f'{"absorbing":>10} {"n_subs":>8} {"cap":>5} {"%abs":>6}')
            if show_subs:
                head += f'  {"med_subs":>9}'
            print(head)
            print('-' * len(head))

            for key in sorted(rows):
                T, m, dT = key
                counts = rows[key]
                n = sum(counts.values())
                pct = 100.0 * counts.get('absorbing', 0) / n if n else 0.0
                flag = '  <-- CHECK' if counts.get('redraw_cap') else ''
                line = (f'{T:>3} {m:>3} {dT:>5.1f} {n:>5} '
                        f'{counts.get("absorbing", 0):>10} '
                        f'{counts.get("n_subs", 0):>8} '
                        f'{counts.get("redraw_cap", 0):>5} {pct:>5.0f}%')
                if show_subs:
                    vals = sorted(subs[key])
                    med = vals[len(vals) // 2] if vals else 0
                    line += f'  {med:>9}'
                print(line + flag)

    total = sum(grand.values())
    if total:
        print(f'\n{"=" * 52}')
        print(f'All conditions shown, {regime}: {total} replicates')
        for reason in ('absorbing', 'n_subs', 'redraw_cap'):
            c = grand.get(reason, 0)
            print(f'  {reason:<12} {c:>8}  {100.0 * c / total:>5.1f}%')


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--cache_dir', default='simulation_cache')
    p.add_argument('--regime', default='simultaneous',
                   choices=['simultaneous', 'sequential', 'all'],
                   help="'simultaneous' is m=T, 'sequential' is m=1.")
    p.add_argument('--dT', type=float, default=None,
                   help='Restrict to one task divergence.')
    p.add_argument('--median_subs', action='store_true',
                   help='Also report the median number of realized '
                        'substitutions.')
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    report(args.cache_dir, args.regime, args.dT, args.median_subs)