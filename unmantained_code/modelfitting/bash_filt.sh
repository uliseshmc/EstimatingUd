#!/bin/bash

# Create logs directory
mkdir -p logs/

# Parameter lists
UPP_LIST=(100 1000 10000)
MAX_LIST=(10 50)
TOL_LIST=(1e-6 1e-8)

# compute sizes
NUM_UPP=${#UPP_LIST[@]}
NUM_MAX=${#MAX_LIST[@]}
NUM_TOL=${#TOL_LIST[@]}
TOTAL=$((NUM_UPP * NUM_MAX * NUM_TOL))

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

#python -m EstimatingUd.modelfitting.testing_trinuc_model -chrm $CHRM -upper $UPP -max_restarts $MAX -tol $TOL
