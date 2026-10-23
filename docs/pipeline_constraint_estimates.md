# Model fitting and constraint

These instructions assume `docs/pipeline_filtering.md` is finished, `Udenv` is active and
`DATA_HUMCHIMPORANGOR114` is set.

## Stage 4: fitting time non-reversible substitution models

```
ud fit --chrm $chr --submodel singlent
ud fit --chrm $chr --submodel trinuc
```

Single-nucleotide fits use cogent3's `GN` model with `time_het="max"` and Orangutan as a
discrete edge. Trinucleotide fits use the models registered by `estimating_ud.models`:
`GNC_CpG_ss` (sense-codon states) for cds with every parameter independent, and
`GT_CpG_ss` (all 64 trinucleotides) for the other regions with `omega` fixed at 1.

The results are pickled to `<region>/chrm<chr>/sm_output/singlent_lh.pickle` or
`trinucleotide_lh.pickle`. These binaries are sensitive to the library versions: fit and
summarise in the same environment.

On the cluster: `scripts/submit.sh scripts/slurm/04_fit.sh`.

## Stage 5: expected number of substitutions and constraint

```
ud constraint
```

For both submodels and chromosomes 1-22 and X, the expected number of substitutions (ENS)
on the Human edge is computed for every region under its own fitted rate matrix and under
the intergenicAR matrix of the same chromosome; constraint is `(ENS_IGAR - ENS)/ENS_IGAR`.
For the trinucleotide cds model the three stop codons, absent from its states, are padded
into the motif probabilities with a 1e-12 pseudocount. intergenicAR appears in the ENS
table as `IGAR`.

Outputs: `results/tables/{singlent,trinuc}_ENS.csv` and `_constraint.csv`. Use
`--submodel`, `--chrm` and `--results-dir` to restrict or redirect; on the cluster,
`scripts/submit.sh scripts/slurm/05_constraint.sh`.

## Summary statistics, mutability and Ud

Run the notebooks in `notebooks/` in order with the `Udenv` kernel:

- `02_summarystats_constraint.ipynb`: mean and SEM of constraint across chromosomes,
  writes `results/tables/*_constraint_mean_sem.csv`.
- `03_estimating_mutrates.ipynb`: applies the Oman et al. 2022 trinucleotide mutability
  model to `human_seq.fa` of every region, writes `results/tables/mutability_motiflength3.csv`.
- `04_summarystats_mutability.ipynb`: normalises to a per-site rate of 1.34e-8, writes
  `average_mutrates_motiflength3_mean_sem.csv`, `UT_perregion_motiflength3.csv` and
  `length_perregion.csv`.
- `05_Ud_estimation.ipynb`: combines the above into Ud and the share of each region.

The `*_mean_sem.csv` files are written by the notebooks with pandas' unstacked layout
(first column `mean`/`sem`, then `Region`, then the value); notebook 05 reads that layout.
