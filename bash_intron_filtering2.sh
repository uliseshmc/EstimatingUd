#!/bin/bash
#SBATCH --job-name=intron_filtering_2
#SBATCH --output=logs/logs_intron_filtering_%A_%a.out
#SBATCH --error=logs/logs_intron_filtering_%A_%a.err
#SBATCH --array=1-24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100gb
#SBATCH --time=3-00:00:00
#SBATCH --account=masel
#SBATCH --partition=standard

# Create logs directory
mkdir -p logs_intron_filtering/

# Initialize conda (using the safer shell.bash hook approach)
eval "$(conda shell.bash hook)"
conda activate UdChimpHumOran

CHROMOSOMES=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"

if (( TASK_ID < 1 || TASK_ID > ${#CHROMOSOMES[@]} )); then
    echo "Invalid SLURM_ARRAY_TASK_ID: $TASK_ID" > &2
    exit 1
fi

CHROMOSOME="${CHROMOSOMES[$((TASK_ID - 1))]}"
echo "Running intron filtering for chromosome $CHROMOSOME"

python3 filtering_gaps_introns_2.py -chrm "$CHROMOSOME"