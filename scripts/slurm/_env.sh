#!/bin/bash
# Sourced by every SLURM wrapper in this folder after the #SBATCH header.
# Activates the single project environment (Udenv), moves to the submission
# directory (the repo root), and maps the array task ID to a chromosome.
set -euo pipefail

module load micromamba/2.0.2-2 2>/dev/null || true
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
MAMBA_EXE="${MAMBA_EXE:-$(command -v micromamba || echo /opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba)}"
eval "$("$MAMBA_EXE" shell hook --shell bash)"
micromamba activate Udenv
command -v ud >/dev/null || { echo "ud is not on PATH: run 'pip install -e .' inside Udenv" >&2; exit 1; }

# sbatch copies the script to its spool dir, so $0 is useless; SLURM_SUBMIT_DIR is where sbatch ran
cd "${SLURM_SUBMIT_DIR:-$PWD}"
: "${DATA_HUMCHIMPORANGOR114:?export DATA_HUMCHIMPORANGOR114=<alignment download root> before sbatch}"

CHROMOSOMES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
TASK_ID="${SLURM_ARRAY_TASK_ID:?must run as an array job (sbatch --array=1-24)}"
if (( TASK_ID < 1 || TASK_ID > ${#CHROMOSOMES[@]} )); then
    echo "Invalid SLURM_ARRAY_TASK_ID: $TASK_ID" >&2
    exit 1
fi
CHROMOSOME="${CHROMOSOMES[$((TASK_ID - 1))]}"
echo "$(date '+%F %T') host=$(hostname) env=Udenv chromosome=$CHROMOSOME data_root=$DATA_HUMCHIMPORANGOR114"
