#!/bin/bash
#SBATCH --job-name=human_seq
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --array=1-24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=14gb
#SBATCH --time=8:00:00
#SBATCH --account=masel
#SBATCH --partition=standard
# Submit from the repo root with scripts/submit.sh (it creates logs/ first).

source "$SLURM_SUBMIT_DIR/scripts/slurm/_env.sh"

ud human-seq --chrm "$CHROMOSOME"
