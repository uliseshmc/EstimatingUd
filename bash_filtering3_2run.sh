#!/bin/bash
#SBATCH --job-name=filtering_3
#SBATCH --output=logs_2run/filtering_3_%A_%a.out
#SBATCH --error=logs_2run/filtering_3_%A_%a.err
#SBATCH --array=0-22
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100gb
#SBATCH --time=3-00:00:00
#SBATCH --account=masel
#SBATCH --partition=standard

# I made this bash because the in the first bash some runs failed
# Recheck runs that failed from first iteration and check that this list is complete
#This list is also available in the file logs/runs_with_errors_rerun.txt

# Create logs directory
mkdir -p logs_2run/

# Initialize conda (using the safer shell.bash hook approach)
eval "$(conda shell.bash hook)"
conda activate UdChimpHumOran

REGIONS=("intergenicAR" "intergenicAR" "intergenicAR" "intergenicAR" "intergenicAR" "intergenicAR" "intergenicAR" "intergenicAR" "intergenicAR" "intergenicAR" "intergenicAR"  "intronsAR" "intronsAR" "intronsAR" "distalIG" "distalIG" "distalIG" "distalIG" "distalIG" "distalIG" "distalIG" "distalIG" "distalIG")
CHROMOSOMES=("1" "2" "3" "4" "5" "6" "8" "9" "11" "12" "X" "4" "12" "X" "4" "5" "13" "14" "15" "6" "17" "18" "19")

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"

REGION="${REGIONS[$TASK_ID]}"
CHROMOSOME="${CHROMOSOMES[$TASK_ID]}"

echo "Running filtering for region=$REGION chromosome=$CHROMOSOME"

python3 filtering_gaps_nointrons_3.py -reg "$REGION" -chrm "$CHROMOSOME"
