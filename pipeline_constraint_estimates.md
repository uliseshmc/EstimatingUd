# Preparation

This instructions assume that pipeline_filtering.md is finished. 
Instructions for installation of the conda package Ensembl0.7.9 can be found there.

Make sure that the variable DATA_HUMCHIMPORANGOR114 in the file paths.py points at the location where you downloaded the alignments. 

Start by activating the environment
```
conda activate Ensembl0.7.9
```

# Fitting time non reversible substitution models

To fit single nucleotide substitution models use 
```
fittinglh_sm_4.py -chrm $chr -submodel singlent
```

To fit trinucleotide substitution models use 
```
fittinglh_sm_4.py -chrm $chr -submodel trinuc
```

These scripts will produce .pickle files inside a /sm_output folder containing binaries of the substitution models. These binaries will be sensitive to changes on the library environment. Make sure to use the same env throughout to avoid mistakes.

# Meassuring contraint estimates

Run
```
estimating_constraint_5.py
```
This script will produce a series of csv files on the output_data/ folder. The files *ENS.csv, *constraint.csv contain raw estimates of the ENS and constraint per genomic region and chromosome. The files *constraint_mean_sem.csv contain the mean and standard deviation accross chromosomes for each region.

# Estimating mutability

Run
```
estimating_mutrates.ipynb
```

# Estimating mutation rates and sequence length

Run
change name!
```
plot_mutability.ipynb
```

# Estimating Uobs

Run
#change name
```
Ud_estimation.ipynb
```

## Bash mode

If using a cluster that runs under a SLURM system, you can use bash*.sh files to run filtering for all regions through all seqids.

Sometimes SLURM fails. Please double checked the .err files to check for common errors. A typical flagged error is AttributeError: 'NotCompleted' object has no attribute 'write'. This usually happens because the length of the alignment is 0. 
To make a list of all other errors I used

```
find . -maxdepth 1 -type f -name "*.err" -size +0c
```
This will output a list of all .err files into runs_with_errors.txt.