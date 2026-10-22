#!/usr/bin/env bash
# First-pass Wright-Fisher pilot. Run from anywhere; uses the repo root.
#   bash wf/run_wf_pilot.sh            # all cores
#   bash wf/run_wf_pilot.sh 32         # 32 workers
# Inside tmux:  tmux new -s wf  ->  bash wf/run_wf_pilot.sh  ->  Ctrl-b d
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p wf_out
WORKERS="${1:-$(nproc)}"
python wf/wf_pilot.py \
    --T 2 4 --dT 1.4 --m 1 T \
    --period 1 10 100 1000 10000 100000 \
    --n_reps 20 --n_gens 2000000 --ckpt 2000 \
    --workers "$WORKERS" 2>&1 | tee -a wf_out/pilot.log
