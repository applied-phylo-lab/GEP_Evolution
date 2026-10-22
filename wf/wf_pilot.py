#!/usr/bin/env python3
"""
wf_pilot.py
===========
Wright-Fisher counterpart of simulate.py's SSWM engine.

Same genotype-phenotype map, performance and fitness as simulate.py (imported,
not re-implemented); same initial genomes and task ensembles per replicate, so
every WF run pairs with its SSWM run in simulation_cache. What changes is the
population genetics:

  * haploid, asexual (fully linked) population of N individuals, tracked as
    genotype classes;
  * every generation: selection + drift by multinomial sampling on fitness
    under the current active task set, then mutation;
  * mutation: each individual mutates with probability mu*L*K per generation,
    flipping one uniformly chosen entry of G (mu per entry, as in SSWM);
    neutral and deleterious mutations are therefore included;
  * the active set of m tasks is redrawn every `period` GENERATIONS, on a fixed
    clock independent of substitutions (in SSWM it was redrawn every
    substitution). For m = T the active set never changes and period is moot.

Recorded at every checkpoint (every `ckpt` generations):
  majority genotype (uint8 L x K), its frequency, number of segregating
  classes, differentiation and optimization of the majority genotype,
  population-mean fitness under the current active set, and the Hamming
  distance of the majority genotype from the initial genome.

Usage (from the repository root):
  python3 wf/wf_pilot.py --bench                       # timing on one replicate
  python3 wf/wf_pilot.py --T 4 --dT 1.4 --m 1 --period 1000 --n_reps 20
  python3 wf/wf_pilot.py --T 2 4 --dT 1.4 --m 1 T --period 1 10 100 1000 10000 100000

Output: wf_out/wf_T{T}_dT{dT}_m{m}_P{period}.npz (+ _meta.json). Conditions
already present with enough replicates are skipped.
"""
import argparse, hashlib, json, os, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed

for _v in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
           'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

import numpy as np
import simulate as S

WF_VERSION = 1


def wf_seed(rep, T, dT, m, period, density):
    key = f'wf_{rep}_{T}_{dT:.10g}_{m}_{period}_{density:.10g}'
    return int(hashlib.md5(key.encode()).hexdigest(), 16) % (2 ** 31)


def task_perf(genome, tasks, gamma):
    """(T,) performances of one genotype on every task."""
    return np.array([S._task_performance(genome, tasks[:, j], gamma)[1]
                     for j in range(tasks.shape[1])])


def class_fitness(Pmat, active, r):
    """Fitness of every class under the active set (uniform weights)."""
    w = np.ones(len(active)) / len(active)
    return S.fitness_power_mean_rows(Pmat[:, active], w, r)


def measure(genome, tasks, gamma, td):
    info = S.compute_performance(genome, tasks, gamma)
    diff = S.mean_pairwise_phenotype_distance(info['phenotype']) / td
    d = info['d']
    opt = 1.0 - float(np.linalg.norm(d)) / np.sqrt(d.shape[0])
    return diff, opt


def run_wf(g0, tasks, T, m, period, n_gens, ckpt, N, mu, gamma, r, seed):
    rng = np.random.default_rng(seed)
    L, K = g0.shape
    td = S.compute_mean_pairwise_distance(tasks)
    U = N * mu * L * K                      # expected new mutations / generation

    genomes = [g0.astype(np.uint8)]
    keys = {g0.astype(np.uint8).tobytes(): 0}
    counts = np.array([N], dtype=np.int64)
    Pmat = task_perf(g0.astype(float), tasks, gamma)[None, :]

    def active_for(epoch):
        if m >= T:
            return np.arange(T)
        return np.sort(rng.choice(T, size=m, replace=False))

    epoch = 0
    active = active_for(epoch)
    W_cur = None                            # fitness of each class; None = stale
    n_ck = n_gens // ckpt + 1
    out = dict(gen=np.zeros(n_ck, np.int64), diff=np.zeros(n_ck),
               opt=np.zeros(n_ck), maj_freq=np.zeros(n_ck),
               n_classes=np.zeros(n_ck, np.int64), mean_W=np.zeros(n_ck),
               hamming0=np.zeros(n_ck, np.int64),
               maj_genome=np.zeros((n_ck, L, K), np.uint8))
    n_switch = 0
    ci = 0

    for g in range(n_gens + 1):
        if g % ckpt == 0:
            j = int(np.argmax(counts))
            G = genomes[j].astype(float)
            out['gen'][ci] = g
            out['diff'][ci], out['opt'][ci] = measure(G, tasks, gamma, td)
            out['maj_freq'][ci] = counts[j] / N
            out['n_classes'][ci] = len(counts)
            W = class_fitness(Pmat, active, r)
            out['mean_W'][ci] = float(np.sum(W * counts) / N)
            out['hamming0'][ci] = int(np.sum(genomes[j] != g0))
            out['maj_genome'][ci] = genomes[j]
            ci += 1
        if g == n_gens:
            break

        if m < T and g > 0 and g % period == 0:
            epoch += 1
            active = active_for(epoch)
            n_switch += 1
            W_cur = None

        # selection + drift (a monomorphic population is unchanged by sampling)
        if len(counts) > 1:
            if W_cur is None:
                W_cur = class_fitness(Pmat, active, r)
            wts = W_cur * counts
            tot = wts.sum()
            p = (wts / tot) if tot > 0 else (counts / counts.sum())
            counts = rng.multinomial(N, p)

        # mutation
        n_mut = rng.poisson(U)
        for _ in range(n_mut):
            par = int(rng.choice(len(counts), p=counts / N))
            child = genomes[par].copy()
            e = int(rng.integers(L * K))
            i, k = divmod(e, K)
            child[i, k] ^= 1
            key = child.tobytes()
            counts[par] -= 1
            if key in keys:
                counts[keys[key]] += 1
            else:
                keys[key] = len(genomes)
                genomes.append(child)
                counts = np.append(counts, 1)
                Pmat = np.vstack([Pmat, task_perf(child.astype(float), tasks, gamma)])
                W_cur = None

        # drop extinct classes
        alive = counts > 0
        if not alive.all():
            idx = np.flatnonzero(alive)
            genomes = [genomes[t] for t in idx]
            counts = counts[idx]
            Pmat = Pmat[idx]
            keys = {gg.tobytes(): t for t, gg in enumerate(genomes)}
            if W_cur is not None:
                W_cur = W_cur[idx]

    out['n_switch'] = n_switch
    out['task_dT_realized'] = td
    return out


def _worker(job):
    (T, dT, m, period, rep, tasks, n_gens, ckpt, N, mu, gamma, r, L, K, dens) = job
    g0 = S.make_genome(L, K, dens, S.make_genome_seed(rep, L, K, dens))
    t0 = time.time()
    res = run_wf(g0, tasks, T, m, period, n_gens, ckpt, N, mu, gamma, r,
                 wf_seed(rep, T, dT, m, period, dens))
    res['seconds'] = time.time() - t0
    return (T, dT, m, period), rep, res


def load_tasks(task_dir, out_dir, L, T, dT, n_reps):
    """Task ensembles identical to the SSWM cache (seeded on rep, T, dT).

        Read from task_dir (the SSWM simulation_cache) when present; otherwise
        built deterministically and stored under out_dir/tasks. Never writes
        into simulation_cache.
    """
    try:
        ens, _, meta = S.load_task_ensembles(S.task_cache_path(task_dir, L, T))
        if dT in ens and meta['n_ensembles'] >= n_reps:
            return ens[dT]
    except (OSError, KeyError, ValueError):
        pass
    own = os.path.join(out_dir, 'tasks')
    os.makedirs(own, exist_ok=True)
    base = os.path.join(own, f'taskens_L{L}_T{T}_dT{dT}')
    try:
        ens, _, meta = S.load_task_ensembles(base)
        if dT in ens and meta['n_ensembles'] >= n_reps:
            return ens[dT]
    except (OSError, KeyError, ValueError):
        pass
    print(f'  building task ensembles T={T} dT={dT} ({n_reps})', flush=True)
    alpha_map, ens, real = S.build_task_ensembles(T, L, [dT], n_reps)
    S.save_task_ensembles(base, ens, alpha_map, real, T, L, n_reps)
    return ens[dT]


def cond_path(out_dir, T, dT, m, period):
    return os.path.join(out_dir, f'wf_T{T}_dT{dT}_m{m}_P{period}')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--T', type=int, nargs='+', default=[4])
    p.add_argument('--dT', type=float, nargs='+', default=[1.4])
    p.add_argument('--m', nargs='+', default=['1', 'T'],
                   help="integers, or 'T' for the simultaneous limit")
    p.add_argument('--period', type=int, nargs='+',
                   default=[1, 10, 100, 1000, 10000, 100000])
    p.add_argument('--n_reps', type=int, default=20)
    p.add_argument('--n_gens', type=int, default=2_000_000)
    p.add_argument('--ckpt', type=int, default=2000)
    p.add_argument('--N', type=int, default=10_000)
    p.add_argument('--mu', type=float, default=1e-7)
    p.add_argument('--L', type=int, default=100)
    p.add_argument('--K', type=int, default=4)
    p.add_argument('--gamma', type=float, default=1.0)
    p.add_argument('--fitness_r', type=float, default=0.0)
    p.add_argument('--density', type=float, default=0.25)
    p.add_argument('--task_dir', default=os.path.join(REPO, 'simulation_cache'),
                   help='where task ensembles are read (or built if missing)')
    p.add_argument('--out_dir', default=os.path.join(REPO, 'wf_out'))
    p.add_argument('--workers', type=int, default=None)
    p.add_argument('--bench', action='store_true',
                   help='time one replicate (T=4, dT=1.4, m=1, period=1000, '
                        '20k generations) and extrapolate')
    a = p.parse_args()

    if a.bench:
        tasks = load_tasks(a.task_dir, a.out_dir, a.L, 4, 1.4, 1)[0]
        g0 = S.make_genome(a.L, a.K, a.density, S.make_genome_seed(0, a.L, a.K, a.density))
        for m in (1, 4):
            t0 = time.time()
            res = run_wf(g0, tasks, 4, m, 1000, 20_000, 1000, a.N, a.mu,
                         a.gamma, a.fitness_r, 1)
            dt = time.time() - t0
            print(f'm={m}: 20k generations in {dt:.1f}s -> '
                  f'{dt / 20_000 * a.n_gens / 60:.1f} min per replicate at '
                  f'n_gens={a.n_gens}; final diff={res["diff"][-1]:.3f}, '
                  f'classes={res["n_classes"][-1]}')
        return

    os.makedirs(a.out_dir, exist_ok=True)
    jobs, todo = [], []
    for T in a.T:
        for dT in a.dT:
            tasks_all = load_tasks(a.task_dir, a.out_dir, a.L, T, dT, a.n_reps)
            ms = sorted({T if x == 'T' else int(x) for x in a.m})
            for m in ms:
                if m > T:
                    continue
                periods = [0] if m == T else a.period
                for P in periods:
                    path = cond_path(a.out_dir, T, dT, m, P)
                    if os.path.exists(path + '_meta.json'):
                        with open(path + '_meta.json') as f:
                            meta = json.load(f)
                        if (meta['n_reps'] >= a.n_reps
                                and meta['n_gens'] >= a.n_gens):
                            print(f'skip (done): {os.path.basename(path)}')
                            continue
                    todo.append((T, dT, m, P))
                    for rep in range(a.n_reps):
                        jobs.append((T, dT, m, max(P, 1), rep, tasks_all[rep],
                                     a.n_gens, a.ckpt, a.N, a.mu, a.gamma,
                                     a.fitness_r, a.L, a.K, a.density))

    print(f'{len(todo)} conditions, {len(jobs)} replicate runs, '
          f'{a.workers or os.cpu_count()} workers', flush=True)
    results = {c: {} for c in todo}
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(_worker, j) for j in jobs]
        for f in as_completed(futs):
            (T, dT, m, P1), rep, res = f.result()
            key = next(c for c in todo if c[:3] == (T, dT, m) and max(c[3], 1) == P1)
            results[key][rep] = res
            if len(results[key]) == a.n_reps:
                T_, dT_, m_, P_ = key
                path = cond_path(a.out_dir, T_, dT_, m_, P_)
                reps = [results[key][i] for i in range(a.n_reps)]
                arrays = {k: np.stack([rr[k] for rr in reps])
                          for k in ('gen', 'diff', 'opt', 'maj_freq',
                                    'n_classes', 'mean_W', 'hamming0',
                                    'maj_genome')}
                arrays['task_dT_realized'] = np.array([rr['task_dT_realized'] for rr in reps])
                arrays['n_switch'] = np.array([rr['n_switch'] for rr in reps])
                arrays['seconds'] = np.array([rr['seconds'] for rr in reps])
                np.savez_compressed(path + '.npz', **arrays)
                with open(path + '_meta.json', 'w') as fh:
                    json.dump(dict(T=T_, dT=dT_, m=m_, period=P_, n_reps=a.n_reps,
                                   n_gens=a.n_gens, ckpt=a.ckpt, N=a.N, mu=a.mu,
                                   L=a.L, K=a.K, gamma=a.gamma,
                                   fitness_r=a.fitness_r, density=a.density,
                                   wf_version=WF_VERSION,
                                   schema_version=S.SCHEMA_VERSION), fh, indent=2)
                d_end = arrays['diff'][:, -1].mean()
                print(f'  saved {os.path.basename(path)}  final diff {d_end:.3f}  '
                      f'({time.time() - t0:.0f}s)', flush=True)
                results.pop(key)
    print(f'Done in {time.time() - t0:.0f}s')


if __name__ == '__main__':
    main()
