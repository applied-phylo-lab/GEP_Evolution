#!/usr/bin/env python3
"""
analyze_wf.py
=============
Summarize wf_out/ against the SSWM cache. numpy only.

For each WF condition: mean differentiation and optimization of the majority
genotype averaged over the last 25% of generations (late window) and at the
halfway point, plus diversity diagnostics. SSWM references for the same
replicates: m = T at its terminal (peak) state; m = 1 at its last state
(substitution 400 in the main cache).

  python3 wf/analyze_wf.py [--wf_dir wf_out] [--cache_dir simulation_cache]
"""
import argparse, glob, json, os
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def sswm_ref(cache_dir, T, dT, m, n):
    pat = os.path.join(cache_dir, 'L100_K4_gamma1.0_fr0.0_v2', 'density0.2500',
                       f'sim_T{T}_dT{dT}_m{m}_alpha*.npz')
    hits = glob.glob(pat)
    if not hits:
        return np.nan
    z = np.load(hits[0])
    meta = json.load(open(hits[0][:-4] + '_meta.json'))['rep_meta']
    return float(np.mean([z[f'rep{i}_pheno_dist'][-1] / meta[i]['task_dT_realized']
                          for i in range(min(n, len(meta)))]))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--wf_dir', default=os.path.join(REPO, 'wf_out'))
    p.add_argument('--cache_dir', default=os.path.join(REPO, 'simulation_cache'))
    a = p.parse_args()

    rows = []
    for mf in sorted(glob.glob(os.path.join(a.wf_dir, 'wf_*_meta.json'))):
        meta = json.load(open(mf))
        z = np.load(mf[:-len('_meta.json')] + '.npz')
        n_ck = z['diff'].shape[1]
        late = slice(int(0.75 * n_ck), n_ck)
        rows.append(dict(
            T=meta['T'], dT=meta['dT'], m=meta['m'], P=meta['period'],
            n=meta['n_reps'],
            diff_half=z['diff'][:, n_ck // 2].mean(),
            diff_late=z['diff'][:, late].mean(),
            diff_late_sd=z['diff'][:, late].mean(axis=1).std(ddof=1),
            opt_late=z['opt'][:, late].mean(),
            maj_freq=z['maj_freq'][:, late].mean(),
            classes=z['n_classes'][:, late].mean(),
            ham=z['hamming0'][:, -1].mean(),
            sswm_mT=sswm_ref(a.cache_dir, meta['T'], meta['dT'], meta['T'], meta['n_reps']),
            sswm_m1=sswm_ref(a.cache_dir, meta['T'], meta['dT'], 1, meta['n_reps']),
            secs=np.median(z['seconds'])))
    rows.sort(key=lambda r: (r['T'], r['dT'], -r['m'] if r['m'] == r['T'] else 0, r['P']))

    print(f'{"T":>2} {"dT":>4} {"m":>2} {"period":>7} {"n":>3} | {"diff@half":>9} '
          f'{"diff late":>9} {"(sd)":>6} {"opt late":>8} | {"maj f":>5} {"classes":>7} '
          f'{"ham0":>5} | {"SSWM m=T":>8} {"SSWM m=1":>8} | {"s/rep":>6}')
    for r in rows:
        P = '-' if r['m'] == r['T'] else str(r['P'])
        print(f'{r["T"]:>2} {r["dT"]:>4} {r["m"]:>2} {P:>7} {r["n"]:>3} | '
              f'{r["diff_half"]:>9.3f} {r["diff_late"]:>9.3f} {r["diff_late_sd"]:>6.3f} '
              f'{r["opt_late"]:>8.3f} | {r["maj_freq"]:>5.2f} {r["classes"]:>7.1f} '
              f'{r["ham"]:>5.0f} | {r["sswm_mT"]:>8.3f} {r["sswm_m1"]:>8.3f} | '
              f'{r["secs"]:>6.0f}')
    print('\ndiff late = mean over the last 25% of generations (majority genotype).'
          '\nSSWM columns: same replicates; m=T at its peak, m=1 at substitution 400'
          ' (still rising for T=4; long-budget plateau was ~0.54 at T=4, ~0.75 at T=2).')


if __name__ == '__main__':
    main()
