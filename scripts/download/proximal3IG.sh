#!/bin/bash
#SBATCH --job-name=proximal3IG
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --array=1-24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --time=2-00:00:00
#SBATCH --account=masel
#SBATCH --partition=standard
# Submit from the repo root with scripts/submit.sh. Needs $DATA_HUMCHIMPORANGOR114/install
# (eti install) and, for the intergenic regions, $DATA_HUMCHIMPORANGOR114/intergenic_coordinates/
# written by notebooks/01_location_inter_intragenic.ipynb.

source "$SLURM_SUBMIT_DIR/scripts/download/_env.sh"

OUT="$DATA/proximal3IG/alldata_chrm${CHRM}"
eti alignments -i "$DATA/install" -od "$OUT" --align_name "10_primates.epo" --ref human \
    --mask "$DATA/allAR_1column.txt" --ref_coords "$DATA/intergenic_coordinates/chrom${CHRM}_proximal3_IG_coordinates.tsv"
mask_report "$OUT"
