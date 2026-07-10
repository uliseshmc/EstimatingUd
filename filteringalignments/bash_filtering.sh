#!/bin/bash
#SBATCH --job-name=eti_alignments
#SBATCH --output=logs/eti_%A_%a.out
#SBATCH --error=logs/eti_%A_%a.err
#SBATCH --array=1-24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20gb                # Explicitly define memory (adjust as needed)
#SBATCH --time=7-00:00:00
#SBATCH --account=masel
#SBATCH --partition=standard

# Create logs directory
mkdir -p logs/

# Initialize conda (using the safer shell.bash hook approach)
eval "$(conda shell.bash hook)"
conda activate UdChimpHumOran

# Array of chromosomes
CHRM_LIST=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)

# Get chromosome name (SLURM_ARRAY_TASK_ID is 1-indexed)
CHRM=${CHRM_LIST[$SLURM_ARRAY_TASK_ID-1]}

echo "Running filtering for Chromosome: $CHRM"

python3 1_splitting_introns.py -chrm ${CHRM} 
python3 2_codon_aligner.py -chrm ${CHRM} 
python3 3_filttering_trinucleotide.py -chrm ${CHRM} 