#!/bin/bash
#SBATCH --job-name=trinuc_filtering_3
#SBATCH --output=logs_2run/trinuc_filtering_3_%A_%a.out
#SBATCH --error=logs_2run/trinuc_filtering_3_%A_%a.err
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

python3 filtering_trinuc_gaps_nointrons_3.py -reg intronsAR -chrm 4
