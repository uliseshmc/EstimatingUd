#!/bin/bash
# Sourced by the eti download wrappers. Same activation as scripts/slurm/_env.sh, plus
# DATA (the alignment root holding install/, intergenic_coordinates/ and the *_1column.txt
# mask lists) and a post-run report on how much of the human sequence was masked.
set -euo pipefail

module load micromamba/2.0.2-2 2>/dev/null || true
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
MAMBA_EXE="${MAMBA_EXE:-$(command -v micromamba || echo /opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba)}"
eval "$("$MAMBA_EXE" shell hook --shell bash)"
micromamba activate Udenv
command -v eti >/dev/null || { echo "eti is not on PATH in Udenv" >&2; exit 1; }

cd "${SLURM_SUBMIT_DIR:-$PWD}"
# fall back to <repo>/.env when the data root is not exported (same rule as estimating_ud.config)
if [ -z "${DATA_HUMCHIMPORANGOR114:-}" ] && [ -f .env ]; then
    set -a; . ./.env; set +a
fi
DATA="${DATA_HUMCHIMPORANGOR114:?export DATA_HUMCHIMPORANGOR114=<alignment download root> before sbatch}"
DATA="${DATA%/}"

CHROMOSOMES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
TASK_ID="${SLURM_ARRAY_TASK_ID:?must run as an array job (sbatch --array=1-24)}"
if (( TASK_ID < 1 || TASK_ID > ${#CHROMOSOMES[@]} )); then
    echo "Invalid SLURM_ARRAY_TASK_ID: $TASK_ID" >&2
    exit 1
fi
CHRM="${CHROMOSOMES[$((TASK_ID - 1))]}"
echo "$(date '+%F %T') host=$(hostname) env=Udenv chromosome=$CHRM data_root=$DATA"
eti --version

# prepare_outdir <outdir>: refuse to extract into a folder that already holds fragments, unless
# the job was submitted with REPLACE_EXISTING=1, in which case eti is told to delete the folder
# first (ETI_FORCE). Mixing two extractions in one folder would silently mix two definitions.
ETI_FORCE=()
prepare_outdir() {
    local dir=$1
    if [ -d "$dir" ] && [ -n "$(find "$dir" -maxdepth 1 -name '*.fa' -print -quit)" ]; then
        if [ "${REPLACE_EXISTING:-0}" != 1 ]; then
            echo "$dir already holds fragments. To DELETE them and extract again, submit with:" >&2
            echo "    REPLACE_EXISTING=1 scripts/submit.sh <this wrapper>" >&2
            exit 1
        fi
        echo "REPLACE_EXISTING=1: eti will delete $dir before writing"
        ETI_FORCE=(--force_overwrite)
    fi
}

# mask_report <outdir> [short]: fraction of fragments with no '?' and the '?' share of the human
# sequence over a sample of fragments. A download whose masking failed shows ~100% files
# without '?' and ~0% masked share (compare: a healthy chromosome has some files without
# '?' because not every fragment overlaps an annotated repeat). Pass "short" for regions made
# of short pieces (utr5, utr3): most of their fragments contain no repeat, so only the masked
# share is informative (chromosome 22: 17% for utr5, 22% for utr3).
mask_report() {
    local dir=$1 n_total n_nomask
    n_total=$(find "$dir" -maxdepth 1 -name '*.fa' | wc -l)
    n_nomask=$(find "$dir" -maxdepth 1 -name '*.fa' -print0 | xargs -0 -r grep -L '?' | wc -l)
    echo "mask report for $dir: $n_nomask of $n_total fragments contain no '?'"
    find "$dir" -maxdepth 1 -name '*.fa' | head -200 | xargs -r awk '
        /^>/ { keep = ($0 ~ /^>homo_sapiens/) ; next }
        keep { masked += gsub(/\?/, ""); total += length($0) }
        END { if (total) printf "mask report: %.1f%% of %d sampled human bases are masked (?)\n", 100*masked/total, total+masked }'
    if [ "${2:-}" = short ]; then
        echo "mask report: short pieces, so many fragments without '?' is expected; judge by the masked share"
    fi
    if [ -d "$dir/not_completed" ]; then
        echo "not_completed records: $(ls "$dir/not_completed" | wc -l)"
    fi
    # always succeed: this is the last command of the wrappers, and a non-zero status here used to
    # mark finished downloads as FAILED in SLURM whenever there was no not_completed folder
    return 0
}
