#!/bin/bash
#SBATCH --job-name=eti_alignments
#SBATCH --output=logs/eti_%A.out
#SBATCH --error=logs/eti_%A.err
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

echo "Running filtering for Chromosome: 22"

python3 1_splitting_introns.py -chrm 22 
python3 2_codon_aligner.py -chrm 22 
python3 3_filttering_trinucleotide.py -chrm 22 