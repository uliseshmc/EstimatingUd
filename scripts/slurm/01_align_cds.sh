#!/bin/bash
#SBATCH --job-name=align_cds
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --array=1-24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --time=4:00:00
#SBATCH --account=masel
#SBATCH --partition=standard
# Submit from the repo root with scripts/submit.sh (it creates logs/ first).

source "$SLURM_SUBMIT_DIR/scripts/slurm/_env.sh"

ud align-cds --chrm "$CHROMOSOME" --submodel singlent
ud align-cds --chrm "$CHROMOSOME" --submodel trinuc
