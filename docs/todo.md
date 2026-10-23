# Open items

Carried over from the 2026-09-24 audit of the script pipeline and the 2026-10-23 restructure review.
Items the restructure resolved are listed at the end so the history stays legible.

## Analysis decisions still open

### The UTR split does not delimit UTRs
In the `introns` alignments both cds and repeats are masked with the same `?` character, in hundreds
of short runs per gene. `sample_UTR5` / `sample_UTR3` / `removeUTRs_fromintrons` (`apps.py`) cut at
the first/last `?`, i.e. at the nearest masked repeat or cds from the transcript ends, typically a
few hundred bp in. Hence introns5UTR and introns3UTR are both ~6.7 Mb genome-wide and nearly equal,
while real 5'UTR exons total ~6-9 Mb and 3'UTR exons ~28-35 Mb. To sample UTRs properly the cds
mask must be distinguishable from the repeat mask (two downloads, or exon coordinates from
`eti dump-genes`).

### Intron alignments with no mask are handled three different ways
When the Human sequence has no `?`: `sample_UTR5` returns the whole alignment, `removeUTRs_fromintrons`
also returns the whole alignment, `sample_UTR3` raises on `None + 1` and is dropped as NotCompleted.
Such an alignment is counted twice. Decide on one behaviour (probably skip it in all three apps).
`tests/test_apps.py` pins the current boundaries so a change is deliberate.

### Two different definitions of "number of sites"
`*_ENS.csv` `aln_length` is the filtered alignment length (three species, ungapped, non-degenerate).
`mutability_motiflength3.csv` `numb_sites` and `length_perregion.csv` are the raw Human sequence
length from `ud human-seq`, with no gap/degenerate filtering and, for cds, no codon alignment or stop
trimming. `05_Ud_estimation.ipynb` multiplies a constraint estimated on the first set by a site count
from the second. Document the choice or unify.

### Mutation-rate means are unweighted across chromosomes
`04_summarystats_mutability.ipynb` normalises per row with site weights, then averages across
chromosomes with a plain mean, so small high-mutability chromosomes are over-weighted. Use a
site-weighted mean across chromosomes or say so.

### The "CpG" scale in both trinucleotide models is bound to omega
`models.py`: `ssym_preds = gn + cpg + omega` then `cpgdecay = ssym_preds[-1]`, which is omega.
`scales["CpG"]` is really "nonsynonymous". Stage 5 uses `expected_number_subs` and never reads the
scales, so current constraint values are unaffected; anything calling `get_scaled_lengths` gets
mislabelled quantities. Fix is `ssym_preds[-2]`, but it is a model change: decide before touching it.

## Environment

- `eti install` with cogent3-h5seqs 0.9.1 (ensembl_tui pins 0.7.3) is unverified. If it fails, see
  the fallback in `docs/pipeline_downloaddata.md`.
- Whether `/opt/ohpc/pub/apps/micromamba` and `~/micromamba/envs/Udenv` are visible on compute nodes
  has not been checked with a job; the wrappers try `module load micromamba/2.0.2-2` first.

## Resolved by the restructure (2026-10-23)

- Repeat masking that failed at download for five regions: re-downloaded 2026-09-24; the download
  wrappers now print a mask report so it cannot go unnoticed again.
- Row accumulation across submodels and swapped alnstat names in stage 5: fixed in
  `stages/constraint.py`; covered by `tests/test_stages.py`.
- Mean/SEM CSV layout: stage 5 no longer writes mean/sem files; the notebooks do, and
  `05_Ud_estimation.ipynb` now filters both `*_mean_sem.csv` files to their mean rows.
- `-mutmotif` flag that was ignored: removed; `ud human-seq` always truncates to whole trinucleotides.
- Stale docs (`DATA_HUMCHIMPORANG115`, `filtered.fa`, `plot_*.ipynb`, coordinate file names),
  wrappers echoing unset variables, the `environment.yml` from an older env, `*.csv` ignored
  everywhere: rewritten.
- The committed trinuc constraint CSVs predated the final chromosome-1 pickles (differences up to
  0.07 for intronsAR chr9); superseded by the full rerun.
