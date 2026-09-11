# GEP_Evolution

Simulation code for *Unicellular and Multicellular Modes of Selection Impose
Distinct Constraints on Cellular Phenotype Evolution* (Kim & Pennell).

A genome is a binary matrix `G` in {0,1}^(L x K) mapping `L` loci onto `K` gene
expression programs. For each of `T` tasks it deploys the non-negative program
combination that best approximates that task's optimum. During each selective
epoch `m` of the `T` task-specific phenotypes contribute jointly to fitness:
`m = 1` is sequential selection, `m = T` simultaneous, and intermediate `m`
interpolates.

## Layout

    simulate.py        model, SSWM engine, task ensembles, cache I/O
    figlib.py          analysis layer: cutoffs, metrics, aggregation, plot grammar
    run_batch.py       run driver; named specs; writes the provenance manifest
    make_data.py       regenerate every simulated condition   (entry point)
    make_figures.py    regenerate every figure                (entry point)
    figures/           one script per manuscript figure
    tests/             52 checks on the engine
    assets/            hand-made figure sources (git-ignored, not reproducible)
    figures_out/       generated figures (git-ignored, reproducible)
    simulation_cache/  trajectories (git-ignored, ~5 GB, rebuild with make_data.py)

Everything in `figures_out/` and `simulation_cache/` can be rebuilt from source.
The hand-made figure sources in `assets/` cannot, and are kept outside version
control; the composited PDFs are distributed with the manuscript.

## Which script makes which figure

| Figure | Script | Command |
|---|---|---|
| 1  | none — hand-made | `assets/F1.pdf` |
| 2  | `figures/fig_tsweep.py`      | `--which F2 --plain_name` |
| 3  | `figures/fig_tsweep.py`      | `--which F3 --plain_name` |
| 4  | `figures/fig_m_by_tasks.py`  | `--plain_name` |
| S1 | `figures/fig_termination.py` | |
| S2 | `figures/fig_regimes.py`     | `--fitness_r -2.0` |
| S3 | `figures/fig_regimes.py`     | `--gamma 4.0` |
| S4 | `figures/compare_K.py`       | `--K 4 6 8 --cutoff 400 --cutoff_scale fixed` |
| S5 | `figures/fig_regimes.py`     | `--density 0.5` |
| S6 | `figures/fig_msweep.py`      | |

Every figure script also accepts `--cutoff_kind exposure`, which compares
conditions after equal selective epochs per task rather than equal
substitutions.

## Reproducing

    python3 run_batch.py --specs baseline --dry_run   # cost estimate
    python3 run_batch.py --specs baseline             # ~5 GB of trajectories
    python3 figures/fig_tsweep.py --which F2 --plain_name

Figure scripts put the repository root on `sys.path` and anchor their default
`--cache_dir` / `--save_dir` at the repository root, so they run from any
working directory.

## Parameter grid

Task numbers are matched on the task-to-program ratio T/K = 0.5, 1, 1.5, 2:

    K = 4   T = 2, 4, 6, 8
    K = 6   T = 3, 6, 9, 12
    K = 8   T = 4, 8, 12, 16

Baseline: L = 100, K = 4, gamma = 1, r = 0, rho = 0.25, N = 1e4, mu = 1e-7,
200 replicates per condition. The K sweep holds rho at 0.25 rather than at 1/K,
so that program number varies without also changing genotype-matrix density.

## License

MIT. See `LICENSE`.

## Requirements

Tested on Python 3.10 with numpy 2.2, scipy 1.15, matplotlib 3.10.
See `requirements.txt`.
