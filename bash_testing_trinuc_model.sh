#!/bin/bash
#SBATCH --job-name=testing_tinuc_model
#SBATCH --output=logstesting_tinuc_model/eti_%A_%a.out
#SBATCH --error=logstesting_tinuc_model/eti_%A_%a.err
#SBATCH --array=1-13
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20gb                # Explicitly define memory (adjust as needed)
#SBATCH --time=7-00:00:00
#SBATCH --account=masel
#SBATCH --partition=standard

# Create logs directory
mkdir -p logstesting_tinuc_model/

# Initialize conda (using the safer shell.bash hook approach)
eval "$(conda shell.bash hook)"
conda activate UdChimpHumOran

# Parameter lists
UPP_LIST=(100 1000 10000)
MAX_LIST=(10 50)
TOL_LIST=(1e-6 1e-8)

# compute sizes
NUM_UPP=${#UPP_LIST[@]}
NUM_MAX=${#MAX_LIST[@]}
NUM_TOL=${#TOL_LIST[@]}
TOTAL=$((NUM_UPP * NUM_MAX * NUM_TOL))

# Validate SLURM_ARRAY_TASK_ID
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
if [ "$TASK_ID" -gt "$TOTAL" ]; then
	echo "Task ID $TASK_ID out of range (total combinations = $TOTAL)" >&2
	exit 1
fi

# zero-based index
IDX=$((TASK_ID - 1))

# compute indices for each list (cartesian product ordering: UPP x MAX x TOL)
BLOCK=$((NUM_MAX * NUM_TOL))
i=$(( IDX / BLOCK ))
rem=$(( IDX % BLOCK ))
j=$(( rem / NUM_TOL ))
k=$(( rem % NUM_TOL ))

UPP=${UPP_LIST[$i]}
MAX=${MAX_LIST[$j]}
TOL=${TOL_LIST[$k]}

CHRM=22
echo "Running filtering for Chromosome: $CHRM; UPP=$UPP MAX=$MAX TOL=$TOL (task $TASK_ID/$TOTAL)"

python -m EstimatingUd.modelfitting.testing_trinuc_model -chrm $CHRM -upper $UPP --maxr $MAX -tol $TOL
