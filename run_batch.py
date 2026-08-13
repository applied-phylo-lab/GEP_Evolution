#!/usr/bin/env python3
"""
run_batch.py
============
Queues the full set of simulation runs behind the manuscript and its
supplementary analyses.

One `simulate.run_simulations` call covers a single (L, K, gamma, fitness_r)
combination, because those four determine the cache root. Genome density is a
subdirectory within that root, so densities sweep inside one call. The specs
below therefore expand to four runs, not five.

  baseline    K=4, gamma=1, r=0, density 1/K, all m
              -> main figures, and the matched-exposure, exposure-diagnostic
                 and collateral-effect supplementary analyses, all of which are
                 re-indexings of these trajectories rather than new simulations
  density     K=4, gamma=1, r=0, density 0.5,  m in {1, T}
  fitness_r   K=4, gamma=1, r=-2, density 1/K, m in {1, T}
  gamma       K=4, gamma=4, r=0,  density 1/K, m in {1, T}
  programs    K=6, gamma=1, r=0,  density 1/K, m in {1, T}

`density` shares a cache root with `baseline` and differs only in the density
subdirectory, so the two are separate calls purely to give them different m
grids. Every robustness spec contrasts sequential against simultaneous
selection and does not need the intermediate m sweep, which dominates cost at
large T; running it anyway would roughly triple those specs.

BLAS threading
--------------
Each worker process is single-threaded work, but NumPy's BLAS will spawn its
own threads inside every worker unless told not to, oversubscribing the node
several-fold. The limits must be set before NumPy is imported, which this
module does at the top. Setting them in the job script as well does no harm.

Runtime estimate
----------------
`calibrate_solve_time` times a single mutant evaluation (genome copy, NNLS
solve, residual, performance) on the machine actually running the batch, and
the planner multiplies it by the work list. The estimate is an UPPER bound on
wall time for two reasons that pull in opposite directions but do not cancel:

  - replicates that reach an absorbing state stop early, often far short of
    their budget, and this is common at m = T;
  - failed epochs force full-repertoire enumeration and are not predictable
    from the work list.

The first dominates in practice, so realized wall time is usually well below
the estimate. Treat it as a ceiling for planning, not a prediction.

Notification
------------
Configure through the environment, never in source:

    export SIM_NOTIFY_EMAIL=you@example.edu
    export SMTP_HOST=smtp.example.edu     # omit to use a local sendmail/mail
    export SMTP_PORT=587
    export SMTP_USER=...                  # omit for unauthenticated relays
    export SMTP_PASSWORD=...
    export SMTP_FROM=...                  # defaults to SIM_NOTIFY_EMAIL

Notification failures are reported and never abort or fail the batch.

Usage:
  python run_batch.py --dry_run
  python run_batch.py
  python run_batch.py --specs baseline gamma
  python run_batch.py --n_reps 100 --workers 64
"""

import os

# Must precede the NumPy import; see "BLAS threading" above.
for _var in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
             'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_var, '1')

import argparse
import copy
import json
import platform
import smtplib
import socket
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone
from email.message import EmailMessage

import numpy as np

import simulate as S


# ============================================================
# 1. SPECS
# ============================================================

BASE_CFG = dict(
    L=100, K=4, GAMMA=1.0, FITNESS_R=0.0,
    T_VALUES=[2, 4, 6, 8],
    M_VALUES=None,
    TASK_DIVERGENCES=[0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4],
    GENOME_DENSITIES=[0.25],
    N_REPS=200, N_POP=10_000, MU=1e-7,
    N_SUBS=None, N_SUBS_EPOCHS=None, N_SUBS_FLOOR=401, N_SUBS_RULE=None,
    MAX_REDRAWS=100, CHUNK_SIZE=10, RECORD_MODULARITY='snapshots',
    TASK_SEED=270, CACHE_DIR='simulation_cache',
    N_WORKERS=None, N_GENOME_SNAPSHOTS=10, TASK_SAMPLING='random',
)

SPECS = {
    'baseline': dict(
        GENOME_DENSITIES=[0.25],
    ),
    'density': dict(
        GENOME_DENSITIES=[0.5],
        M_VALUES='endpoints',
    ),
    'fitness_r': dict(
        FITNESS_R=-2.0,
        M_VALUES='endpoints',
    ),
    'gamma': dict(
        GAMMA=4.0,
        M_VALUES='endpoints',
    ),
    'programs': dict(
        K=6,
        M_VALUES='endpoints',
    ),
}

SPEC_ORDER = ['baseline', 'density', 'fitness_r', 'gamma', 'programs']

SPEC_NOTES = {
    'baseline': 'main figures, and every re-indexed supplementary analysis',
    'density': 'denser initial genotype matrix',
    'fitness_r': 'negative power mean, isolating fitness aggregation',
    'gamma': 'steeper performance decline, isolating curvature',
    'programs': 'six programs, testing whether the T<=K boundary moves with K',
}


def build_cfg(spec_name: str, overrides: dict) -> dict:
    cfg = copy.deepcopy(BASE_CFG)
    cfg.update(SPECS[spec_name])
    cfg.update({k: v for k, v in overrides.items() if v is not None})
    return cfg


# ============================================================
# 2. ENVIRONMENT
# ============================================================

def git_revision() -> str:
    try:
        out = subprocess.run(['git', 'rev-parse', 'HEAD'],
                             capture_output=True, text=True, timeout=5)
        if out.returncode == 0:
            rev = out.stdout.strip()
            dirty = subprocess.run(['git', 'status', '--porcelain'],
                                   capture_output=True, text=True, timeout=5)
            return rev + ('-dirty' if dirty.stdout.strip() else '')
    except (OSError, subprocess.SubprocessError):
        pass
    return 'unknown'


def environment_record() -> dict:
    return {
        'host': socket.gethostname(),
        'python': sys.version.split()[0],
        'numpy': np.__version__,
        'platform': platform.platform(),
        'git': git_revision(),
        'schema_version': S.SCHEMA_VERSION,
        'cpu_count': os.cpu_count(),
        'blas_threads': {v: os.environ.get(v) for v in
                         ('OMP_NUM_THREADS', 'MKL_NUM_THREADS',
                          'OPENBLAS_NUM_THREADS')},
    }


def default_workers() -> int:
    """Usable worker count. Prefers the CPU affinity mask over os.cpu_count(),
    since the mask reflects any restriction placed on the process (cgroup,
    taskset, container) while the count reports the whole machine."""
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


def calibrate_solve_time(L: int, K: int, n: int = 2000, seed: int = 0) -> float:
    """Seconds per mutant evaluation on this machine.

    The unit counted by `simulate.estimate_cost` is one mutant-task evaluation:
    copy the genome, flip a bit, solve the NNLS, take the residual and the
    performance. Timing that whole operation rather than the bare solver keeps
    the estimate honest about copy and norm overhead, which is not negligible
    at these matrix sizes.
    """
    rng = np.random.default_rng(seed)
    genome = (rng.random((L, K)) < 0.25).astype(float)
    task = rng.dirichlet(np.ones(L))
    task /= np.linalg.norm(task)

    for _ in range(50):                                    # warm up
        S._task_performance(genome, task, 1.0)

    t0 = time.perf_counter()
    for i in range(n):
        mut = genome.copy()
        mut[i % L, i % K] = 1.0 - mut[i % L, i % K]
        S._task_performance(mut, task, 1.0)
    return (time.perf_counter() - t0) / n


def format_hours(seconds: float) -> str:
    if seconds < 3600:
        return f'{seconds / 60:.0f} min'
    if seconds < 48 * 3600:
        return f'{seconds / 3600:.1f} h'
    return f'{seconds / 86400:.1f} d'


# ============================================================
# 3. NOTIFICATION
# ============================================================

def send_email(subject: str, body: str, to_addr: str) -> str:
    """Send via SMTP if SMTP_HOST is set, otherwise hand off to a local
    sendmail or mail binary. Returns a short status string; never raises."""
    if not to_addr:
        return 'skipped (no SIM_NOTIFY_EMAIL)'

    host = os.environ.get('SMTP_HOST')
    from_addr = os.environ.get('SMTP_FROM', to_addr)

    try:
        if host:
            msg = EmailMessage()
            msg['Subject'] = subject
            msg['From'] = from_addr
            msg['To'] = to_addr
            msg.set_content(body)

            port = int(os.environ.get('SMTP_PORT', '587'))
            user = os.environ.get('SMTP_USER')
            password = os.environ.get('SMTP_PASSWORD')

            with smtplib.SMTP(host, port, timeout=30) as smtp:
                if port != 25:
                    smtp.starttls()
                if user and password:
                    smtp.login(user, password)
                smtp.send_message(msg)
            return f'sent via {host}'

        for binary in (['sendmail', '-t'], ['mail', '-s', subject, to_addr]):
            try:
                payload = (f'To: {to_addr}\nSubject: {subject}\n\n{body}'
                           if binary[0] == 'sendmail' else body)
                subprocess.run(binary, input=payload, text=True, check=True,
                               capture_output=True, timeout=60)
                return f'sent via {binary[0]}'
            except FileNotFoundError:
                continue
        return 'failed (no SMTP_HOST, and no sendmail or mail binary)'

    except Exception as exc:                       # notification is best-effort
        return f'failed ({type(exc).__name__}: {exc})'


# ============================================================
# 4. BATCH
# ============================================================

def cost_of(cfg: dict) -> dict:
    """Plan one spec without running it. Loads or builds task ensembles, since
    the work list cannot be expanded without them."""
    alpha_maps, ensembles = {}, {}
    for T in cfg['T_VALUES']:
        a, e, _ = S.ensure_task_ensembles(
            cfg['CACHE_DIR'], cfg['L'], T,
            cfg['TASK_DIVERGENCES'], cfg['N_REPS'])
        alpha_maps[T], ensembles[T] = a, e

    work = S.build_work(cfg, alpha_maps, ensembles)
    total, longest = S.estimate_cost(work)
    return {
        'work_items': len(work),
        'conditions': len({(w[0], w[1], w[2], w[3]) for w in work}),
        'total_solves': total,
        'longest_item_solves': longest,
        'budgets': sorted({(w[0], w[2], w[7]) for w in work}),
    }


def run_batch(spec_names, overrides, dry_run, manifest_path, notify_to):
    started = datetime.now(timezone.utc)
    t0 = time.time()
    manifest = {
        'started_utc': started.isoformat(),
        'environment': environment_record(),
        'specs': {},
    }

    n_workers = overrides.get('N_WORKERS') or default_workers()
    print(f'Batch of {len(spec_names)} spec(s): {", ".join(spec_names)}')
    print(f'Workers: {n_workers}   '
          f'BLAS threads: {os.environ.get("OMP_NUM_THREADS")}')

    solve_times = {}
    for K in sorted({build_cfg(n, overrides)['K'] for n in spec_names}):
        L = BASE_CFG['L']
        solve_times[K] = calibrate_solve_time(L, K)
        print(f'Calibration: L={L} K={K}  '
              f'{solve_times[K] * 1e6:.1f} us per mutant evaluation')
    manifest['solve_time_us'] = {k: v * 1e6 for k, v in solve_times.items()}
    manifest['n_workers'] = n_workers
    print()

    grand_total = 0.0
    grand_eta = 0.0
    for name in spec_names:
        cfg = build_cfg(name, overrides)
        print(f'--- {name}: {SPEC_NOTES.get(name, "")}')
        print(f'    K={cfg["K"]} gamma={cfg["GAMMA"]} r={cfg["FITNESS_R"]} '
              f'densities={cfg["GENOME_DENSITIES"]} m={cfg["M_VALUES"]}')

        entry = {'cfg': {k: v for k, v in cfg.items() if k != 'N_SUBS_RULE'}}
        try:
            plan = cost_of(cfg)
            per_solve = solve_times[cfg['K']]
            eta = plan['total_solves'] * per_solve / n_workers
            floor = plan['longest_item_solves'] * per_solve
            plan['eta_seconds'] = eta
            plan['critical_path_seconds'] = floor
            entry['plan'] = {k: v for k, v in plan.items() if k != 'budgets'}
            grand_total += plan['total_solves']
            grand_eta += eta
            print(f'    {plan["conditions"]} conditions, '
                  f'{plan["total_solves"] / 1e9:.2f}e9 solves')
            print(f'    ETA <= {format_hours(eta)}   '
                  f'(critical path {format_hours(floor)})')
            if dry_run:
                for T, m, n_subs in plan['budgets']:
                    print(f'      T={T:>2d} m={m:>2d} n_subs={n_subs}')
        except Exception:
            entry['status'] = 'plan_failed'
            entry['traceback'] = traceback.format_exc()
            print(f'    PLAN FAILED\n{entry["traceback"]}')
            manifest['specs'][name] = entry
            continue

        if dry_run:
            entry['status'] = 'planned'
            manifest['specs'][name] = entry
            continue

        spec_t0 = time.time()
        try:
            S.run_simulations(cfg)
            entry['status'] = 'ok'
        except Exception:
            entry['status'] = 'failed'
            entry['traceback'] = traceback.format_exc()
            print(f'    SPEC FAILED\n{entry["traceback"]}')
        entry['wall_seconds'] = round(time.time() - spec_t0, 1)
        manifest['specs'][name] = entry
        print(f'    {entry["status"]} in {entry["wall_seconds"]:.0f}s\n')

    manifest['finished_utc'] = datetime.now(timezone.utc).isoformat()
    manifest['wall_seconds'] = round(time.time() - t0, 1)
    manifest['total_solves'] = grand_total
    manifest['eta_seconds'] = grand_eta
    manifest['dry_run'] = dry_run

    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2, default=str)
    print(f'Manifest: {manifest_path}')

    statuses = {k: v.get('status') for k, v in manifest['specs'].items()}
    if dry_run:
        print(f'\nBatch ETA <= {format_hours(grand_eta)} on {n_workers} workers')
        print('Upper bound: replicates reaching an absorbing state stop early, '
              'which is common at m=T.')
    else:
        print(f'Total wall time: {format_hours(manifest["wall_seconds"])} '
              f'(estimated <= {format_hours(grand_eta)})')
    print(f'Status: {statuses}')

    if notify_to and not dry_run:
        ok = all(v == 'ok' for v in statuses.values())
        subject = (f'[sim] batch {"complete" if ok else "FINISHED WITH ERRORS"} '
                   f'on {manifest["environment"]["host"]}')
        lines = [
            f'Started  {manifest["started_utc"]}',
            f'Finished {manifest["finished_utc"]}',
            f'Wall     {manifest["wall_seconds"] / 3600:.2f} h',
            f'Host     {manifest["environment"]["host"]}',
            f'Git      {manifest["environment"]["git"]}',
            '',
            'Specs:',
        ]
        for name, entry in manifest['specs'].items():
            lines.append(f'  {name:<12} {entry.get("status")}'
                         f'  {entry.get("wall_seconds", "-")}s')
            if entry.get('traceback'):
                lines.append('    ' + entry['traceback'].strip().splitlines()[-1])
        lines += ['', f'Manifest: {os.path.abspath(manifest_path)}']
        body = '\n'.join(lines)
        print(f'Notification: {send_email(subject, body, notify_to)}')

    return manifest


# ============================================================
# 5. CLI
# ============================================================

def parse_args():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[2])
    p.add_argument('--specs', nargs='+', default=SPEC_ORDER,
                   choices=SPEC_ORDER,
                   help='Which specs to run, in the order given.')
    p.add_argument('--dry_run', action='store_true',
                   help='Plan and cost every spec without running anything.')
    p.add_argument('--n_reps', type=int, default=None)
    p.add_argument('--workers', type=int, default=None)
    p.add_argument('--chunk_size', type=int, default=None)
    p.add_argument('--cache_dir', type=str, default=None)
    p.add_argument('--manifest', type=str, default=None,
                   help='Manifest path; defaults to a timestamped file.')
    p.add_argument('--email', type=str, default=None,
                   help='Override SIM_NOTIFY_EMAIL.')
    p.add_argument('--no_email', action='store_true')
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()

    overrides = {
        'N_REPS': args.n_reps,
        'N_WORKERS': args.workers or default_workers(),
        'CHUNK_SIZE': args.chunk_size,
        'CACHE_DIR': args.cache_dir,
    }

    manifest_path = args.manifest or (
        f'batch_manifest_{datetime.now().strftime("%Y%m%d_%H%M%S")}.json')

    notify_to = None if args.no_email else (
        args.email or os.environ.get('SIM_NOTIFY_EMAIL'))

    run_batch(args.specs, overrides, args.dry_run, manifest_path, notify_to)