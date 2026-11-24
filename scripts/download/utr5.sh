#!/bin/bash
#SBATCH --job-name=utr5
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --array=1-24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --time=06:00:00
#SBATCH --account=masel
#SBATCH --partition=standard
# Submit from the repo root with scripts/submit.sh. Needs $DATA_HUMCHIMPORANGOR114/install
# (eti install) and $DATA_HUMCHIMPORANGOR114/genic_coordinates/ written by `ud genic-coords`.
# 5' UTR exons of each protein-coding gene's canonical transcript, all repeats masked.
# Chromosome 22 took about 15 minutes (2026-11-24); chromosome 1 has about five times as many pieces.

source "$SLURM_SUBMIT_DIR/scripts/download/_env.sh"

OUT="$DATA/utr5/alldata_chrm${CHRM}"
prepare_outdir "$OUT"
eti alignments -i "$DATA/install" -od "$OUT" --align_name "10_primates.epo" --ref human \
    --mask "$DATA/allAR_1column.txt" --ref_coords "$DATA/genic_coordinates/chrom${CHRM}_utr5.tsv" \
    ${ETI_FORCE[@]+"${ETI_FORCE[@]}"}
mask_report "$OUT" short
