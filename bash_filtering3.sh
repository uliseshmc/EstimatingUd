#!/bin/bash
#SBATCH --job-name=filtering_3
#SBATCH --output=logs/filtering_3_%A_%a.out
#SBATCH --error=logs/filtering_3_%A_%a.err
#SBATCH --array=1-120
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100gb
#SBATCH --time=3-00:00:00
#SBATCH --account=masel
#SBATCH --partition=standard

# Create logs directory
mkdir -p logs/

# Initialize conda (using the safer shell.bash hook approach)
eval "$(conda shell.bash hook)"
conda activate UdChimpHumOran

REGIONS=("intergenicAR" "intronsAR" "distalIG" "proximal5IG" "proximal3IG")
CHROMOSOMES=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
TOTAL_COMBINATIONS=$(( ${#REGIONS[@]} * ${#CHROMOSOMES[@]} ))

if (( TASK_ID < 1 || TASK_ID > TOTAL_COMBINATIONS )); then
    echo "Invalid SLURM_ARRAY_TASK_ID: $TASK_ID" >&2
    exit 1
fi

REGION_INDEX=$(((TASK_ID - 1) / ${#CHROMOSOMES[@]}))
CHROMOSOME_INDEX=$(((TASK_ID - 1) % ${#CHROMOSOMES[@]}))
REGION="${REGIONS[$REGION_INDEX]}"
CHROMOSOME="${CHROMOSOMES[$CHROMOSOME_INDEX]}"

echo "Running filtering for region=$REGION chromosome=$CHROMOSOME"

python3 filtering_gaps_nointrons_3.py -reg "$REGION" -chrm "$CHROMOSOME"
