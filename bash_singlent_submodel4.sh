#!/bin/bash
#SBATCH --job-name=singlent_submodel_4
#SBATCH --output=logs_submodel/singlent_submodel_4_%A_%a.out
#SBATCH --error=logs_submodel/singlent_submodel_4_%A_%a.err
#SBATCH --array=1-216
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=200gb
#SBATCH --time=2-00:00:00
#SBATCH --account=masel
#SBATCH --partition=standard

# Create logs directory
mkdir -p logs_submodel/

# Initialize conda (using the safer shell.bash hook approach)
eval "$(conda shell.bash hook)"
conda activate Ensembl0.7.9

CHROMOSOMES=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"

if (( TASK_ID < 1 || TASK_ID > ${#CHROMOSOMES[@]} )); then
    echo "Invalid SLURM_ARRAY_TASK_ID: $TASK_ID" >&2
    exit 1
fi

CHROMOSOME="${CHROMOSOMES[$((TASK_ID - 1))]}"

echo "Running substitution model for region=$REGION chromosome=$CHROMOSOME"

python3 fittinglh_singlentmodel_4.py -chrm "$CHROMOSOME" -submodel singlent
python3 fittinglh_singlentmodel_4.py -chrm "$CHROMOSOME" -submodel trinuc