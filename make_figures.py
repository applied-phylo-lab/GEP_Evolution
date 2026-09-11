#!/usr/bin/env python3
"""
make_figures.py
===============
Regenerate every figure in the manuscript and the supplement, under the names
the LaTeX source expects.

Each figure script prints `Saved: <path>`; this runner reads that line and moves
the result to its canonical name, so figures_out holds only F1-F4 and FS1-FS6
and the figure-to-file mapping lives in one place. Running a figure script
directly still produces its parameterised name, which is the safe default
because it cannot overwrite another parameter set. Figure 1 is hand-made and is
copied from assets/.

Requires the simulation cache; build it with `python3 make_data.py`.

Usage:
  python3 make_figures.py
  python3 make_figures.py --only F2 F3 F4
  python3 make_figures.py --list
"""
import argparse, io, os, re, shutil, subprocess, sys, time

REPO = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(REPO, 'figures_out')

# canonical name -> (script or None, arguments, one-line description)
FIGURES = [
    ('F1',  None,                  [],
     'schematic (hand-made, copied from assets/)'),
    ('F2',  'fig_tsweep.py',       ['--which', 'F2'],
     'sequential selection, m = 1, two cutoffs'),
    ('F3',  'fig_tsweep.py',       ['--which', 'F3'],
     'simultaneous selection, m = T, two cutoffs'),
    ('F4',  'fig_m_by_tasks.py',   [],
     'm = 1, T/2, T against task number, one column per dT'),
    ('FS1', 'fig_termination.py',  [],
     'stationarity and absorbing states'),
    ('FS2', 'fig_regimes.py',      ['--fitness_r', '-2.0'],
     'negative power-mean fitness, r = -2'),
    ('FS3', 'fig_regimes.py',      ['--gamma', '4.0'],
     'sharper performance function, gamma = 4'),
    ('FS4', 'compare_K.py',        ['--K', '4', '6', '8', '--cutoff', '400',
                                    '--cutoff_scale', 'fixed'],
     'program number, T/K on the x-axis'),
    ('FS5', 'fig_regimes.py',      ['--density', '0.5'],
     'denser initialisation, rho = 0.5'),
    ('FS6', 'fig_msweep.py',       [],
     'gain over sequential selection against m'),
]
SAVED = re.compile(r'^Saved:\s*(.+\.\w+)\s*$', re.M)


def run_one(name, script, extra, verbose):
    if script is None:
        src = os.path.join(REPO, 'assets', 'F1.pdf')
        dst = os.path.join(OUT, 'F1.pdf')
        shutil.copy2(src, dst)
        return dst
    path = os.path.join(REPO, 'figures', script)
    # Not every script defines --no_summary; ask the script itself rather than
    # keeping a list here that will drift.
    src = io.open(path, encoding='utf-8').read()
    quiet = ['--no_show'] + (['--no_summary'] if "'--no_summary'" in src else [])
    cmd = [sys.executable, path] + quiet + extra
    r = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write(r.stdout[-3000:] + r.stderr[-3000:])
        raise SystemExit(f'{name}: {script} exited {r.returncode}')
    if verbose:
        print(r.stdout.rstrip()[-2000:])
    hits = SAVED.findall(r.stdout)
    if not hits:
        raise SystemExit(f'{name}: {script} printed no "Saved:" line')
    produced = hits[-1]
    canonical = os.path.join(OUT, f'{name}.pdf')
    if os.path.abspath(produced) != os.path.abspath(canonical):
        # move, so figures_out only ever holds F1-F4 and FS1-FS6. Running a
        # figure script directly still gives the parameterised name, which is
        # the safe default because it cannot overwrite another parameter set.
        shutil.move(produced, canonical)
    return canonical


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--only', nargs='+', default=None)
    p.add_argument('--list', action='store_true')
    p.add_argument('--verbose', action='store_true')
    args = p.parse_args()

    if args.list:
        for name, script, extra, desc in FIGURES:
            print(f'  {name:<4} {script or "(hand-made)":<20} {desc}')
        return

    os.makedirs(OUT, exist_ok=True)
    todo = [f for f in FIGURES if args.only is None or f[0] in args.only]
    t0 = time.time()
    for name, script, extra, desc in todo:
        t = time.time()
        path = run_one(name, script, extra, args.verbose)
        print(f'{name:<4} {os.path.basename(path):<12} '
              f'{time.time() - t:6.1f}s   {desc}', flush=True)
    print(f'\n{len(todo)} figures in {time.time() - t0:.0f}s -> {OUT}')


if __name__ == '__main__':
    main()
