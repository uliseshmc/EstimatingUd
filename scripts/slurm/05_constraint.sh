#!/bin/bash
#SBATCH --job-name=constraint
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=30gb
#SBATCH --time=4:00:00
#SBATCH --account=masel
#SBATCH --partition=standard
# Submit from the repo root with scripts/submit.sh. Not an array job: one run covers both
# submodels and chromosomes 1-22 + X, writing results/tables/{singlent,trinuc}_{ENS,constraint}.csv
set -euo pipefail
export SLURM_ARRAY_TASK_ID=1   # satisfies _env.sh; the chromosome it picks is unused here
source "$SLURM_SUBMIT_DIR/scripts/slurm/_env.sh"
ud constraint
