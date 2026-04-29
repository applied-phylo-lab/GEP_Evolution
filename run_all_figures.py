#!/usr/bin/env python3
"""
run_all_figures.py
==================
Master script that runs all figure scripts in sequence.

Generates:
  fig_2               — pheno_dist + usage at fixed cutoff, m=1 vs m=T, x=T/K
  fig_3               — pheno_dist vs m sweep, one panel per T/K, at fixed cutoff
  figS_endpoint_metrics — 3-row endpoint summary (pheno, cosine, Jaccard) vs m
  figS_trajectories   — full substitution trajectories of dZ/dT
  figS_tradeoff_modularity — mutational tradeoff-modularity plane trajectories

All figures are saved to --save_dir with gamma and fitness_r in the filename
so runs with different parameters never overwrite each other.

Usage:
  python run_all_figures.py
  python run_all_figures.py --gamma 2.0 --fitness_r -1.0 --cutoff 100
  python run_all_figures.py --densities 0.25 0.5 --no_show
"""

import argparse
import os

import matplotlib.pyplot as plt

import fig_2
import fig_3
import figS_endpoint_metrics
import figS_trajectories
import figS_tradeoff_modularity


# ============================================================
# CLI
# ============================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Generate all figures from simulation cache."
    )
    # Shared simulation parameters
    p.add_argument("--cache_dir",  type=str,   default="simulation_cache")
    p.add_argument("--save_dir",   type=str,   default="figures_out")
    p.add_argument("--fmt",        type=str,   default="pdf")
    p.add_argument("--L",          type=int,   default=100)
    p.add_argument("--K",          type=int,   default=4)
    p.add_argument("--gamma",      type=float, default=1.0)
    p.add_argument("--fitness_r",  type=float, default=0.0,
                   help="Power mean exponent: 0=geometric mean, >0 soft-max, <0 soft-min")
    p.add_argument("--densities",  type=float, nargs="+", default=[0.25])
    p.add_argument("--T",          type=int,   nargs="+",
                   default=[2, 3, 4, 6, 8], dest="T_values")
    p.add_argument("--dT",         type=float, nargs="+",
                   default=[0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4], dest="task_divs")
    # Cutoff applies to fig_2, fig_3, figS_endpoint_metrics
    p.add_argument("--cutoff",     type=int,   default=200,
                   help="Substitution step cutoff for endpoint figures (default 200)")
    # figS_tradeoff_modularity extras
    p.add_argument("--plane_cache_dir", type=str, default="plane_cache")
    p.add_argument("--force_recompute", action="store_true",
                   help="Force recompute of tradeoff-modularity plane cache")
    p.add_argument("--no_show",    action="store_true")
    p.add_argument("--no_summary", action="store_true")
    return p.parse_args()


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    args = parse_args()
    os.makedirs(args.save_dir, exist_ok=True)

    shared = dict(
        cache_dir=args.cache_dir,
        L=args.L,
        K=args.K,
        gamma=args.gamma,
        fitness_r=args.fitness_r,
        T_values=args.T_values,
        task_divs=args.task_divs,
    )

    # ----------------------------------------------------------
    # fig_2
    # ----------------------------------------------------------
    print("\n" + "=" * 60)
    print("fig_2")
    print("=" * 60)
    cfg2 = fig_2.DataConfig(
        **shared,
        densities=args.densities,
        cutoffs=[args.cutoff],
    )
    out2 = fig_2.OutputConfig(
        save_dir=args.save_dir, filename="fig_2", fmt=args.fmt,
        show=not args.no_show, print_summary=not args.no_summary,
    )
    fig_cfg2 = fig_2.FigureConfig()

    for density in cfg2.densities:
        data, task_maps, _ = fig_2.load_all_for_density(cfg2, density)
        for cutoff in cfg2.cutoffs:
            def sp2(tag):
                return os.path.join(
                    out2.save_dir,
                    f"{out2.filename}_{tag}_sub{cutoff}"
                    f"_gamma{cfg2.gamma}_fr{cfg2.fitness_r}"
                    f"_density{density:.4f}.{out2.fmt}"
                )
            f_main = fig_2.fig_main_plot(data, task_maps, cfg2, fig_cfg2,
                                          cutoff=cutoff, save_path=sp2("main"))
            f_usage = fig_2.fig_usage_subset(data, cfg2, fig_cfg2,
                                              cutoff=cutoff, save_path=sp2("usage"))
            if not args.no_summary:
                fig_2.print_summary(data, task_maps, cfg2, fig_cfg2, cutoff)
            if args.no_show:
                plt.close(f_main)
                plt.close(f_usage)
            else:
                plt.show()

    # ----------------------------------------------------------
    # fig_3
    # ----------------------------------------------------------
    print("\n" + "=" * 60)
    print("fig_3")
    print("=" * 60)
    cfg3 = fig_3.DataConfig(
        **shared,
        densities=args.densities,
        cutoff=args.cutoff,
    )
    out3 = fig_3.OutputConfig(
        save_dir=args.save_dir, filename="fig_3", fmt=args.fmt,
        show=not args.no_show, print_summary=not args.no_summary,
    )
    fig_cfg3 = fig_3.FigureConfig()

    for density in cfg3.densities:
        data, task_maps, _ = fig_3.load_all_for_density(cfg3, density)
        cutoff_tag = f"_sub{cfg3.cutoff}" if cfg3.cutoff is not None else ""
        save_path = os.path.join(
            out3.save_dir,
            f"{out3.filename}{cutoff_tag}"
            f"_gamma{cfg3.gamma}_fr{cfg3.fitness_r}"
            f"_density{density:.4f}.{out3.fmt}"
        )
        f3 = fig_3.fig3_m_sweep(data, task_maps, cfg3, fig_cfg3, save_path=save_path)
        if not args.no_summary:
            fig_3.print_summary(data, task_maps, cfg3)
        if args.no_show:
            plt.close(f3)
        else:
            plt.show()

    # ----------------------------------------------------------
    # figS_endpoint_metrics
    # ----------------------------------------------------------
    print("\n" + "=" * 60)
    print("figS_endpoint_metrics")
    print("=" * 60)
    cfgS = figS_endpoint_metrics.DataConfig(
        **shared,
        densities=args.densities,
        cutoff=args.cutoff,
    )
    outS = figS_endpoint_metrics.OutputConfig(
        save_dir=args.save_dir, filename="figS_endpoint_metrics", fmt=args.fmt,
        show=not args.no_show, print_summary=not args.no_summary,
    )
    fig_cfgS = figS_endpoint_metrics.FigureConfig()

    for density in cfgS.densities:
        data, task_maps, _ = figS_endpoint_metrics.load_all_for_density(cfgS, density)
        cutoff_tag = f"_sub{cfgS.cutoff}" if cfgS.cutoff is not None else ""
        save_path = os.path.join(
            outS.save_dir,
            f"{outS.filename}{cutoff_tag}"
            f"_gamma{cfgS.gamma}_fr{cfgS.fitness_r}"
            f"_density{density:.4f}.{outS.fmt}"
        )
        fS = figS_endpoint_metrics.fig2_m_sweep(
            data=data, task_maps=task_maps,
            data_cfg=cfgS, fig_cfg=fig_cfgS, save_path=save_path,
        )
        if not args.no_summary:
            figS_endpoint_metrics.print_summary(data, task_maps, cfgS)
        if args.no_show:
            plt.close(fS)
        else:
            plt.show()

    # ----------------------------------------------------------
    # figS_trajectories
    # ----------------------------------------------------------
    print("\n" + "=" * 60)
    print("figS_trajectories")
    print("=" * 60)
    cfgT = figS_trajectories.DataConfig(
        **shared,
        densities=args.densities,
    )
    outT = figS_trajectories.OutputConfig(
        save_dir=args.save_dir, filename="figS_trajectories", fmt=args.fmt,
        show=not args.no_show,
    )
    fig_cfgT = figS_trajectories.FigureConfig()

    for density in cfgT.densities:
        data, task_maps, _ = figS_trajectories.load_all_for_density(cfgT, density)
        tag = f"_gamma{cfgT.gamma}_fr{cfgT.fitness_r}_density{density:.4f}"

        fT1 = figS_trajectories.fig_trajectories(
            data=data, task_maps=task_maps,
            data_cfg=cfgT, fig_cfg=fig_cfgT,
            save_path=os.path.join(outT.save_dir,
                                   f"{outT.filename}{tag}.{outT.fmt}"),
        )
        fT2a = figS_trajectories.fig_trajectories_m_sweep(
            data=data, task_maps=task_maps,
            data_cfg=cfgT, fig_cfg=fig_cfgT, target_tk=0.8,
            save_path=os.path.join(outT.save_dir,
                                   f"{outT.filename}_m_sweep_tk0.8{tag}.{outT.fmt}"),
        )
        fT2b = figS_trajectories.fig_trajectories_m_sweep(
            data=data, task_maps=task_maps,
            data_cfg=cfgT, fig_cfg=fig_cfgT, target_tk=2.0,
            save_path=os.path.join(outT.save_dir,
                                   f"{outT.filename}_m_sweep_tk2.0{tag}.{outT.fmt}"),
        )
        if args.no_show:
            for f_ in [fT1, fT2a, fT2b]:
                if f_ is not None:
                    plt.close(f_)
        else:
            plt.show()

    # ----------------------------------------------------------
    # figS_tradeoff_modularity  (single density per run)
    # ----------------------------------------------------------
    print("\n" + "=" * 60)
    print("figS_tradeoff_modularity")
    print("=" * 60)
    for density in args.densities:
        cfgM = figS_tradeoff_modularity.DataConfig(
            **shared,
            plane_cache_dir=args.plane_cache_dir,
            density=density,
        )
        os.makedirs(args.plane_cache_dir, exist_ok=True)

        save_path = os.path.join(
            args.save_dir,
            f"figS_tradeoff_modularity"
            f"_gamma{cfgM.gamma}_fr{cfgM.fitness_r}"
            f"_density{density:.4f}.{args.fmt}"
        )
        fM = figS_tradeoff_modularity.make_figure(
            data_cfg=cfgM,
            fig_cfg=figS_tradeoff_modularity.FigureConfig(),
            save_path=save_path,
            force_recompute=args.force_recompute,
        )
        if args.no_show:
            plt.close(fM)
        else:
            plt.show()

    print("\nAll figures done.")
