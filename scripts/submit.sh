#!/bin/bash
# Usage, from the repo root:  scripts/submit.sh scripts/slurm/04_fit.sh [extra sbatch args]
# Creates logs/ (SLURM opens the --output file before the job script runs, so the
# directory must exist at submission time), checks the data root is set, then sbatch.
set -euo pipefail
script=${1:?usage: scripts/submit.sh <wrapper.sh> [sbatch args]}; shift
[ -f pyproject.toml ] || { echo "run from the repo root" >&2; exit 1; }
: "${DATA_HUMCHIMPORANGOR114:?export DATA_HUMCHIMPORANGOR114=<alignment download root>}"
mkdir -p logs
sbatch --export=ALL "$@" "$script"
