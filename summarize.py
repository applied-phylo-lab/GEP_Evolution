#!/usr/bin/env python3
"""
summarize.py
============
Generates numerical summaries of simulation results as CSV tables.

Tables produced:
  summary_endpoint.csv     : mean ± std of endpoint metrics per (T, dT, m)
  summary_trajectory.csv   : mean eff_rank at 10%, 25%, 50%, 75%, 100% of subs
  summary_performance.csv  : mean ± std of per-task performance at endpoint
  summary_pheno_dist.csv   : mean ± std of pheno_dist at endpoint

Run:
  python summarize.py
  python summarize.py --cache_dir simulation_cache --L 100 --K 6
"""

import argparse
import csv
import json
import os
from typing import Dict, List, Optional

import numpy as np

from simulate import load_condition, load_task_map, sim_cache_path, task_cache_path


# ============================================================
# DEFAULTS
# ============================================================

DEFAULTS = dict(
    CACHE_DIR='simulation_cache',
    L=100,
    K=6,
    GAMMA=1.0,
    DENSITY=None,
    T_VALUES=[3, 6, 9],
    TASK_DIVERGENCES=[0.2, 0.6, 1.2],
    SAVE_DIR='summary_out',
)


# ============================================================
# DATA LOADING  (same as figures.py)
# ============================================================

def load_all_conditions(
    cache_dir: str, L: int, K: int, gamma: float,
    density: float, T_values: List[int],
    task_divergences: List[float],
) -> Dict:
    data = {}
    for T in T_values:
        data[T] = {}
        tpath = task_cache_path(cache_dir, L, K, gamma, T)
        tpath_meta = tpath + '_meta.json'
        if not os.path.exists(tpath_meta):
            print(f'  WARNING: task map meta not found for T={T}')
            continue
        with open(tpath_meta) as f:
            meta = json.load(f)
        alpha_map = {float(k): v for k, v in meta['alpha_map'].items()}

        for dT in task_divergences:
            data[T][dT] = {}
            if dT not in alpha_map:
                continue
            alpha = alpha_map[dT]

            m = 1
            while m <= T:
                spath = sim_cache_path(
                    cache_dir, L, K, gamma, density, T, dT, m, alpha)
                if os.path.exists(spath + '.npz'):
                    results, _ = load_condition(spath)
                    data[T][dT][m] = results
                m += 1

        for dT in task_divergences:
            if data[T][dT]:
                ms = sorted(data[T][dT].keys())
                print(f'  T={T} dT={dT}: m = {ms}')

    return data


# ============================================================
# SUMMARY FUNCTIONS
# ============================================================

def endpoint_summary(data: Dict, T_values: List[int],
                     task_divergences: List[float]) -> List[dict]:
    """
    For each (T, dT, m): mean ± std of eff_rank, pheno_dist,
    mean_performance, modularity_entropy at final substitution step.
    """
    rows = []
    for T in T_values:
        for dT in task_divergences:
            if dT not in data[T]:
                continue
            for m in sorted(data[T][dT].keys()):
                reps = data[T][dT][m]
                eff_ranks, pheno_dists, mean_perfs, modularities = [], [], [], []
                n_actual = []

                for r in reps:
                    n = r['n_actual_subs']
                    n_actual.append(n)
                    eff_ranks.append(float(r['eff_rank'][n - 1]))
                    pheno_dists.append(float(r['pheno_dist'][n - 1]))
                    mean_perfs.append(float(np.mean(r['P'][n - 1])))
                    mod = r['modularity_entropy'][n - 1]
                    modularities.append(float(mod) if not np.isnan(mod) else np.nan)

                rows.append({
                    'T': T, 'dT': dT, 'm': m,
                    'n_reps': len(reps),
                    'mean_n_subs': float(np.mean(n_actual)),
                    'eff_rank_mean': float(np.mean(eff_ranks)),
                    'eff_rank_std':  float(np.std(eff_ranks)),
                    'pheno_dist_mean': float(np.mean(pheno_dists)),
                    'pheno_dist_std':  float(np.std(pheno_dists)),
                    'mean_perf_mean': float(np.mean(mean_perfs)),
                    'mean_perf_std':  float(np.std(mean_perfs)),
                    'modularity_mean': float(np.nanmean(modularities)),
                    'modularity_std':  float(np.nanstd(modularities)),
                })
    return rows


def trajectory_summary(data: Dict, T_values: List[int],
                       task_divergences: List[float],
                       quantiles: List[float] = None) -> List[dict]:
    """
    For each (T, dT, m): mean eff_rank at fractional substitution
    steps (e.g., 10%, 25%, 50%, 75%, 100% of n_subs).
    """
    if quantiles is None:
        quantiles = [0.10, 0.25, 0.50, 0.75, 1.00]

    rows = []
    for T in T_values:
        for dT in task_divergences:
            if dT not in data[T]:
                continue
            for m in sorted(data[T][dT].keys()):
                reps = data[T][dT][m]
                row = {'T': T, 'dT': dT, 'm': m}

                for q in quantiles:
                    vals = []
                    for r in reps:
                        n = r['n_actual_subs']
                        idx = max(0, min(int(q * n) - 1, n - 1))
                        vals.append(float(r['eff_rank'][idx]))
                    col = f'eff_rank_at_{int(q*100)}pct'
                    row[col + '_mean'] = float(np.mean(vals))
                    row[col + '_std']  = float(np.std(vals))

                rows.append(row)
    return rows


def per_task_performance_summary(data: Dict, T_values: List[int],
                                  task_divergences: List[float]) -> List[dict]:
    """
    For each (T, dT, m): mean ± std of each individual task's
    performance at endpoint. Shows whether performance is balanced
    across tasks or concentrated on some.
    """
    rows = []
    for T in T_values:
        for dT in task_divergences:
            if dT not in data[T]:
                continue
            for m in sorted(data[T][dT].keys()):
                reps = data[T][dT][m]
                # collect (n_reps, T) array of endpoint performances
                P_end = []
                for r in reps:
                    n = r['n_actual_subs']
                    P_end.append(r['P'][n - 1])
                P_arr = np.array(P_end)  # (n_reps, T)

                row = {'T': T, 'dT': dT, 'm': m}
                for t in range(T):
                    row[f'P_task{t+1}_mean'] = float(P_arr[:, t].mean())
                    row[f'P_task{t+1}_std']  = float(P_arr[:, t].std())

                # also report std across tasks (imbalance measure)
                row['P_std_across_tasks_mean'] = float(
                    np.mean([np.std(P_end[i]) for i in range(len(P_end))]))

                rows.append(row)
    return rows


def m_comparison_summary(data: Dict, T_values: List[int],
                          task_divergences: List[float]) -> List[dict]:
    """
    Direct comparison of m=1 vs m=T endpoint metrics.
    Computes difference and effect size (Cohen's d) for eff_rank.
    """
    rows = []
    for T in T_values:
        for dT in task_divergences:
            if dT not in data[T]:
                continue
            m_vals = sorted(data[T][dT].keys())
            m1  = min(m_vals)
            mT  = max(m_vals)

            if m1 not in data[T][dT] or mT not in data[T][dT]:
                continue

            def endpoint_arr(m, key):
                reps = data[T][dT][m]
                return np.array([float(r[key][r['n_actual_subs'] - 1])
                                 for r in reps])

            er_m1 = endpoint_arr(m1, 'eff_rank')
            er_mT = endpoint_arr(mT, 'eff_rank')
            pd_m1 = endpoint_arr(m1, 'pheno_dist')
            pd_mT = endpoint_arr(mT, 'pheno_dist')

            # Cohen's d
            pooled_std = np.sqrt((er_m1.std()**2 + er_mT.std()**2) / 2)
            cohens_d = ((er_mT.mean() - er_m1.mean()) / pooled_std
                        if pooled_std > 1e-12 else np.nan)

            rows.append({
                'T': T, 'dT': dT,
                'm1': m1, 'mT': mT,
                'eff_rank_m1_mean': er_m1.mean(),
                'eff_rank_m1_std':  er_m1.std(),
                'eff_rank_mT_mean': er_mT.mean(),
                'eff_rank_mT_std':  er_mT.std(),
                'eff_rank_diff':    er_mT.mean() - er_m1.mean(),
                'eff_rank_cohens_d': cohens_d,
                'pheno_dist_m1_mean': pd_m1.mean(),
                'pheno_dist_mT_mean': pd_mT.mean(),
                'pheno_dist_diff':    pd_mT.mean() - pd_m1.mean(),
            })
    return rows


# ============================================================
# CSV WRITER
# ============================================================

def write_csv(rows: List[dict], path: str):
    if not rows:
        print(f'  No data for {path}')
        return
    fieldnames = list(rows[0].keys())
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (f'{v:.6f}' if isinstance(v, float) else v)
                             for k, v in row.items()})
    print(f'  Saved: {path}  ({len(rows)} rows)')


# ============================================================
# MAIN
# ============================================================

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--cache_dir', type=str, default=DEFAULTS['CACHE_DIR'])
    p.add_argument('--L',         type=int, default=DEFAULTS['L'])
    p.add_argument('--K',         type=int, default=DEFAULTS['K'])
    p.add_argument('--gamma',     type=float, default=DEFAULTS['GAMMA'])
    p.add_argument('--density',   type=float, default=DEFAULTS['DENSITY'])
    p.add_argument('--T',         type=int, nargs='+',
                   default=DEFAULTS['T_VALUES'], dest='T_VALUES')
    p.add_argument('--dT',        type=float, nargs='+',
                   default=DEFAULTS['TASK_DIVERGENCES'],
                   dest='TASK_DIVERGENCES')
    p.add_argument('--save_dir',  type=str, default=DEFAULTS['SAVE_DIR'])
    return p.parse_args()


if __name__ == '__main__':
    args   = parse_args()
    L      = args.L
    K      = args.K
    gamma  = args.gamma
    density = args.density if args.density is not None else 1.0 / K
    T_values  = args.T_VALUES
    task_divs = args.TASK_DIVERGENCES
    save_dir  = args.save_dir

    os.makedirs(save_dir, exist_ok=True)

    print('Loading data...')
    data = load_all_conditions(
        args.cache_dir, L, K, gamma, density,
        T_values, task_divs
    )

    print('\nGenerating summaries...')

    write_csv(
        endpoint_summary(data, T_values, task_divs),
        os.path.join(save_dir, 'summary_endpoint.csv')
    )

    write_csv(
        trajectory_summary(data, T_values, task_divs),
        os.path.join(save_dir, 'summary_trajectory.csv')
    )

    write_csv(
        per_task_performance_summary(data, T_values, task_divs),
        os.path.join(save_dir, 'summary_per_task_performance.csv')
    )

    write_csv(
        m_comparison_summary(data, T_values, task_divs),
        os.path.join(save_dir, 'summary_m_comparison.csv')
    )

    print(f'\nDone. Tables saved to {save_dir}/')
