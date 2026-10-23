#!/bin/bash
#SBATCH --job-name=filtering_quality
#SBATCH --output=logs_filtering_quality/filtering_quality_%A_%a.out
#SBATCH --error=logs_filtering_quality/filtering_quality_%A_%a.err
#SBATCH --array=1-24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=60gb
#SBATCH --time=4:00:00
#SBATCH --account=masel
#SBATCH --partition=standard

# Create logs directory
mkdir -p logs_filtering_quality/

# Initialize conda (using the safer shell.bash hook approach)
eval "$(conda shell.bash hook)"
conda activate cog2026.10.22_En0.7.9

CHROMOSOMES=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"

if (( TASK_ID < 1 || TASK_ID > ${#CHROMOSOMES[@]} )); then
    echo "Invalid SLURM_ARRAY_TASK_ID: $TASK_ID" >&2
    exit 1
fi

CHROMOSOME="${CHROMOSOMES[$((TASK_ID - 1))]}"

echo "Running alignment-quality filtering for chromosome=$CHROMOSOME"

python3 filtering_alignment_quality_IGAR.py -chrm "$CHROMOSOME"