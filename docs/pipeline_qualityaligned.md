# Alignment quality: fitting the best-aligned windows

These instructions assume `docs/pipeline_downloaddata.md` is finished (`Udenv` active,
`DATA_HUMCHIMPORANGOR114` set, raw `eti` output under `<region>/alldata_chrm<chr>/`). The
branch does not depend on stages 1 to 5; it reads the raw intergenicAR fragments directly.

The main pipeline (`docs/pipeline_filtering.md`) removes masked and gapped columns but keeps
every remaining column whatever the local alignment quality. Mis-aligned stretches look like
excess divergence and inflate the neutral rate matrix that every constraint is measured
against. This branch keeps only the windows of intergenicAR whose sum-of-pairs score beats a
threshold and refits the trinucleotide model on them, so the fit can be compared with the
stage 4 intergenicAR fit. It is exploratory: stage 5 never reads its output.

Three `ud quality` commands, implemented in `src/estimating_ud/stages/quality.py`, run in
order per chromosome; a fourth summarises the fits of every chromosome in one table. A
fifth, `masked-concat`, is a side product described at the end.

```
ud quality masked-databank --chrm $chr
ud quality filter-spscore  --chrm $chr --window 60 --threshold 0.5
ud quality fit-spscore     --chrm $chr --window 60 --threshold 0.5
ud quality ens-spscore                 --window 60 --threshold 0.5
```

`--window` and `--threshold` default to `regions.QUALITY_WINDOW` (60) and
`regions.QUALITY_THRESHOLD` (0.5).

## Step 1: one gap-preserving file per fragment

```
ud quality masked-databank --chrm $chr
```

For intergenicAR, intronsAR and distalIG (`regions.MASKED_REGIONS`), every fragment in
`<region>/alldata_chrm<chr>/` is renamed to Human/Chimpanzee/Orangutan, then

- `omit_degenerates(gap_is_degen=False, motif_length=1)` removes columns holding a masked
  `?` or a non-ACGT symbol but keeps gaps, which `sp_score` needs for its gap penalties;
- `omit_gap_pos(allowed_frac=0.90)` removes the columns that are gaps in every species. With
  three sequences a column is 0, 1/3, 2/3 or 3/3 gaps, so 0.90 only catches the all-gap
  columns; on chromosome 22 these are 59-69% of the raw columns and `sp_score` skips them
  anyway, so dropping them makes every window the same number of scorable columns.

Each filtered fragment is written to `<region>/nonmasked_alldata_chrm<chr>/` with its
original file name. Fragments that end up empty go to the data store's `not_completed`
list instead of a file. The raw `alldata_chrm<chr>` folders are read only.

## Step 2: keep the windows above the sp_score threshold

```
ud quality filter-spscore --chrm $chr --window 60 --threshold 0.5
```

Only intergenicAR (`regions.QUALITY_REGIONS`) is processed. Each file in
`intergenicAR/nonmasked_alldata_chrm<chr>/` goes through `apps.filter_windows_spscore` and
then `omit_degenerates(motif_length=3)`; the completed results are concatenated into

```
intergenicAR/quality_aligned_chrm<chr>/trinucleotide_filtered_window<W>_qualitythreshold<T>.fa
```

`filter_windows_spscore` cuts the fragment into non-overlapping windows of `--window`
columns and scores each with cogent3's `sp_score` (`calc="JC69"`, `gap_insert=12`,
`gap_extend=1`). The score turns each pairwise distance into an estimated number of
unchanged positions, subtracts the affine gap penalty, and sums over the three species pairs,
so three identical gap-free sequences score `3 * window` and gap-dominated windows go
negative. The threshold is not given as a score but as a fraction, `--threshold`: it is the
`sp_score` of a synthetic alignment with the same species in which two sequences are all `A`
and the third has `round(window * threshold)` bases followed by gaps. With the defaults the
threshold is 36 against a maximum of 180. Windows scoring above it are kept and joined end
to end; the trailing partial window is always dropped; windows `sp_score` cannot score are
dropped; a fragment with no passing window returns `NotCompleted` and is skipped.

`omit_degenerates(motif_length=3)` then removes every trinucleotide with a gap, so the file
is comparable to the `trinucleotide_filtered.fa` of stage 3. Keep `--window` a multiple of 3:
the kept windows are concatenated before the trinucleotide filter, and any other width
shifts the reading frame at every window boundary.

No `_alnstat.txt` is written for this file; its length is `len(aln)` after loading it.

## Step 3: fit the trinucleotide model

```
ud quality fit-spscore --chrm $chr --window 60 --threshold 0.5
```

`--window` and `--threshold` only select the input file from step 2. The fit is the stage 4
non-cds fit (`stages.fit.trinuc_model_noncds`): `GT_CpG_ss` with every parameter independent
across edges, `omega` fixed at 1 and motif probabilities not optimised, so the result is
directly comparable with `intergenicAR/chrm<chr>/sm_output/trinucleotide_lh.pickle`. The
result is pickled to

```
intergenicAR/quality_aligned_chrm<chr>/sm_output/trinucleotide_lh_window<W>_qualitythreshold<T>.pickle
```

so fits with different settings sit side by side. As in stage 4, the pickle must be read in
the environment that wrote it.

## Step 4: ENS of every fit, and the comparison with the main pipeline

```
ud quality ens-spscore --window 60 --threshold 0.5
```

Runs once, not per chromosome. For chromosomes 1-22 and X it loads the step 3 pickle and,
with the helpers stage 5 uses (`stages.constraint.human_edge`, cogent3's
`expected_number_subs` at `t=1`), writes the Human-edge ENS to

```
results/tables/trinuc_ENS_spscore_window<W>_qualitythreshold<T>.csv
```

with the stage 5 ENS columns (`Region` is `IGAR`, `Chromosome` is text) plus `window` and
`threshold`. `aln_length` is `len(result.alignment)` of the pickled result, since step 2
writes no `_alnstat.txt`. A missing pickle is an error (`--chrm` restricts the set); Y is
never included because no fragment of it keeps a window (see SLURM below). Stage 5 does not
read this table.

`notebooks/exploratory/ens_quality_aligned.ipynb` merges it with the `IGAR` rows of
`results/tables/trinuc_ENS.csv` (the stage 4 fit of every window, same model), tabulates
the per-chromosome ratio of the two ENS values and of the two alignment lengths, and plots
the mean +/- SEM across chromosomes of both ENS values as two bars
(`results/figures/ens_igar_quality_vs_all_window<W>_qualitythreshold<T>.png`). The
notebook's first cell says what the comparison is: the model is identical, but the stage 3
filter (whole trinucleotides dropped around every `?` or gap in the raw frame) and the
step 1-2 filters (`?` columns dropped singly, all-gap columns removed, re-framed, windows
selected, then gapped trinucleotides dropped) keep different column sets, and the default
threshold essentially removes indel-dense windows, not mismatch-dense ones.

The exploratory notebooks in `notebooks/exploratory/` support the choice of settings:
`plotting_alignment_quality.ipynb` plots the per-window `sp_score` distribution of the
chromosome 22 step 1 files for several window sizes (figures in
`results/figures/alignment_quality_spscore_*.png`), and `app_alignment_quality.ipynb` checks
that `filter_windows_spscore` reproduces the clean-window subset selected there. The
defaults, 60 columns and 0.5, come from that chromosome 22 exploration.

## SLURM

One `--array=1-24` wrapper per step, submitted from the repo root in order, waiting for each
array to finish before the next:

```
scripts/submit.sh scripts/slurm/q1_masked_databank.sh
scripts/submit.sh scripts/slurm/q2_filter_spscore.sh
scripts/submit.sh scripts/slurm/q3_fit_spscore.sh
```

The wrappers hard-code `--window 60 --threshold 0.5`; edit `q2` and `q3` together when
changing them. Memory and time requested: 8 GB / 4 h, 4 GB / 4 h and 40 GB / 8 h (the
2026-10-26 run peaked below 3 GB and 4 min per chromosome for the fit). Logs go to
`logs/<job name>_<array job>_<task>.out|err`; the `NotCompleted` error described in
`docs/pipeline_filtering.md` has the same meaning here (no fragment of the chromosome kept
any window).

## Side product: `masked-concat`

```
ud quality masked-concat --chrm $chr
```

Applies the step 1 filter to the same three regions but concatenates the fragments into
`<region>/chrm<chr>/singlent_nomasked.fa` with a `singlent_nomasked_alnstat.txt`. Nothing in
the pipeline or the notebooks reads it; it was used to inspect the gap content of whole
chromosomes in the archived chromosome 22 analysis. There is no SLURM wrapper for it.
