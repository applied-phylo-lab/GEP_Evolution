# Wright-Fisher pilot

Tests whether the SSWM contrast between sequential (m = 1) and simultaneous
(m = T) selection survives full population-genetic dynamics, and how it depends
on how fast the selected phenotype switches.

`wf_pilot.py` reuses simulate.py's genotype-phenotype map, performance and
fitness, and the same initial genomes and task ensembles per replicate, so each
WF run pairs with its SSWM run. Differences from SSWM:

* haploid asexual population of N individuals, tracked as genotype classes;
* selection + drift by multinomial sampling every generation, then mutation
  (mu per entry of G per generation; neutral and deleterious mutations included);
* for m < T, the active task set is redrawn every `period` generations on a
  fixed clock. In SSWM it was redrawn after every substitution, so the SSWM
  regime corresponds roughly to a period near the substitution interval
  (~10^2-10^3 generations here).

Expectation at the fast end: when the period is much shorter than a sweep,
long-run lineage growth is the product of fitness across environments, i.e. the
geometric mean over tasks -- the m = T fitness. m = 1 should converge to m = T.

## Running (HPC)

    conda activate /local/storage/mk2687/GEP_Evolution/env
    python wf/wf_pilot.py --bench          # timing; then
    tmux new -s wf
    bash wf/run_wf_pilot.sh                # Ctrl-b d to detach
    python wf/analyze_wf.py                # summary table

Task ensembles are read from simulation_cache/tasks_L100_v2/ when present
(copy taskens_T2* and taskens_T4* there for exact pairing); otherwise they are
rebuilt deterministically into wf_out/tasks/ (identical to ~1e-16).

Output: wf_out/wf_T{T}_dT{dT}_m{m}_P{period}.npz + _meta.json, one per
condition. Finished conditions are skipped on re-run.
