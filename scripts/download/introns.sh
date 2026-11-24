#!/bin/bash
#SBATCH --job-name=introns
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --array=1-24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --time=1-00:00:00
#SBATCH --account=masel
#SBATCH --partition=standard
# Submit from the repo root with scripts/submit.sh. Needs $DATA_HUMCHIMPORANGOR114/install
# (eti install) and $DATA_HUMCHIMPORANGOR114/genic_coordinates/ written by `ud genic-coords`.
# Introns: the merged protein-coding gene spans with every canonical-transcript exon (mrna, cds)
# and all repeats masked. exons_allAR_1column.txt is `mrna` followed by cds_allAR_1column.txt:
#     { echo mrna; cat cds_allAR_1column.txt; } > exons_allAR_1column.txt
# Until 2026-11 this wrapper extracted one fragment per gene with only cds and repeats masked;
# a folder still holding those fragments must be replaced, not added to (see prepare_outdir).
# 40 gene regions of chromosome 22 took 8 minutes (2026-11-24); chromosome 1 is estimated at 8-9 hours.

source "$SLURM_SUBMIT_DIR/scripts/download/_env.sh"

OUT="$DATA/introns/alldata_chrm${CHRM}"
MASK="$DATA/exons_allAR_1column.txt"
if [ ! -f "$MASK" ]; then
    echo "missing $MASK; build it with: { echo mrna; cat cds_allAR_1column.txt; } > exons_allAR_1column.txt" >&2
    exit 1
fi
prepare_outdir "$OUT"
eti alignments -i "$DATA/install" -od "$OUT" --align_name "10_primates.epo" --ref human \
    --mask "$MASK" --ref_coords "$DATA/genic_coordinates/chrom${CHRM}_genes_merged.tsv" \
    ${ETI_FORCE[@]+"${ETI_FORCE[@]}"}
mask_report "$OUT"
