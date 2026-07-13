#!/bin/bash
#SBATCH --job-name=filtering_3
#SBATCH --output=logs_3run/filtering_3_%A_%a.out
#SBATCH --error=logs_3run/filtering_3_%A_%a.err
#SBATCH --array=0-6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=140gb
#SBATCH --time=3-00:00:00
#SBATCH --account=masel
#SBATCH --partition=standard

# I made this bash because the in the second bash some runs failed
# Recheck runs that failed from second iteration and check that this list is complete
#This list is also available in the file logs_2run/runs_with_errors_rerun.txt

# Create logs directory
mkdir -p logs_3run/

# Initialize conda (using the safer shell.bash hook approach)
eval "$(conda shell.bash hook)"
conda activate UdChimpHumOran

REGIONS=("intergenicAR" "intronsAR" "intronsAR" "intronsAR" "distalIG" "distalIG" "distalIG")
CHROMOSOMES=("3" "1" "5" "6" "6" "15" "16")

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"

REGION="${REGIONS[$TASK_ID]}"
CHROMOSOME="${CHROMOSOMES[$TASK_ID]}"

echo "Running filtering for region=$REGION chromosome=$CHROMOSOME"

python3 filtering_gaps_nointrons_3.py -reg "$REGION" -chrm "$CHROMOSOME"
