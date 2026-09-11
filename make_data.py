#!/usr/bin/env python3
"""
make_data.py
============
Generate every simulation condition used by the manuscript and the supplement.

This is the reproducibility entry point. `run_batch.py`'s SPECS cover the
baseline and three robustness axes, but the two program-number sweeps were
originally launched with command-line grid overrides that SPECS does not
record, so a clone of the repository could not reproduce them. They are written
out explicitly below.

Task numbers are matched on the task-to-program ratio T/K = 0.5, 1, 1.5, 2:

    K = 4   T = 2, 4, 6, 8      K = 6   T = 3, 6, 9, 12
    K = 8   T = 4, 8, 12, 16

Conditions already cached at sufficient depth are skipped by run_batch.py, so
this is safe to re-run and cheap to resume.

Usage:
  python3 make_data.py --dry_run      # cost estimate for the whole paper
  python3 make_data.py                # run everything (days on one machine)
  python3 make_data.py --only baseline programs_K8
"""
import argparse, os, subprocess, sys, time

REPO = os.path.dirname(os.path.abspath(__file__))

# name -> arguments appended to `python3 run_batch.py`
JOBS = {
    'baseline':    ['--specs', 'baseline'],
    'density':     ['--specs', 'density'],
    'fitness_r':   ['--specs', 'fitness_r'],
    'gamma':       ['--specs', 'gamma'],
    'programs_K6': ['--specs', 'programs', '--K', '6',
                    '--T', '3', '6', '9', '12'],
    'programs_K8': ['--specs', 'programs', '--K', '8',
                    '--T', '4', '8', '12', '16'],
}
ORDER = ['baseline', 'density', 'fitness_r', 'gamma',
         'programs_K6', 'programs_K8']

FIGURES_OF = {
    'baseline':    'Figures 2, 3, 4; S1; S6',
    'density':     'S5',
    'fitness_r':   'S2',
    'gamma':       'S3',
    'programs_K6': 'S4',
    'programs_K8': 'S4',
}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--only', nargs='+', default=ORDER, choices=ORDER)
    p.add_argument('--dry_run', action='store_true')
    p.add_argument('--workers', type=int, default=None)
    p.add_argument('--n_reps', type=int, default=None)
    args = p.parse_args()

    t0 = time.time()
    for name in [j for j in ORDER if j in args.only]:
        cmd = [sys.executable, os.path.join(REPO, 'run_batch.py')] + JOBS[name]
        if args.dry_run:
            cmd.append('--dry_run')
        if args.workers:
            cmd += ['--workers', str(args.workers)]
        if args.n_reps:
            cmd += ['--n_reps', str(args.n_reps)]
        print(f'\n{"=" * 72}\n{name}   -> {FIGURES_OF[name]}\n'
              f'{" ".join(cmd[1:])}\n{"=" * 72}', flush=True)
        r = subprocess.run(cmd, cwd=REPO)
        if r.returncode != 0:
            sys.exit(f'{name} failed with exit code {r.returncode}')
    print(f'\nAll requested data complete in {time.time() - t0:.0f}s.')


if __name__ == '__main__':
    main()
