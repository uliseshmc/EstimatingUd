# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

A research pipeline (not a distributable library) for estimating Ud, the deleterious mutation rate,
from Human/Chimpanzee/Orangutan whole-genome alignments (Ensembl 114, downloaded with `ensembl_tui`,
CLI `eti`) using `cogent3` substitution models. The analysis code is an installable package with a
Click CLI (`ud`); it runs per chromosome, usually as SLURM array jobs, on multi-terabyte alignment
data stored outside the repo. Read the three documents in `docs/` before changing download,
filtering or model-fitting logic.

## Layout

- `src/estimating_ud/` is the only code that matters.
  - `apps.py`: cogent3 composable apps (species renamers, UTR samplers, Human-sequence extractors,
    motif counter, `filter_windows_spscore`).
  - `models.py`: `GNC_CpG_ss` (codon states, cds) and `GT_CpG_ss` (64 trinucleotide states, non-cds),
    registered with cogent3 when the package is imported. `dinuc_models.py` holds the `GDN_CpG*`
    dinucleotide models. **Keep the two modules separate**: both define `make_gn_preds` and friends
    with different bodies, and merging them silently changes the trinucleotide models.
  - `regions.py`: chromosome / submodel / region lists and every on-disk path helper.
  - `config.py`: `data_root()` from `DATA_HUMCHIMPORANGOR114`, falling back to `<repo>/.env`
    (gitignored; `.env.example` is the template) so notebooks work however the kernel was started;
    `results_dir()` resolved from the repo (override with `ESTIMATING_UD_RESULTS`), never from the CWD.
  - `stages/`: one module per stage with a plain `run(...)` function; `cli.py` only parses options.
- `scripts/slurm/` (stage wrappers) and `scripts/download/` (eti wrappers) are `--array=1-24` jobs;
  `scripts/submit.sh` creates `logs/` and calls sbatch. Both `_env.sh` files activate `Udenv`,
  `cd "$SLURM_SUBMIT_DIR"`, and refuse to run without `SLURM_ARRAY_TASK_ID` and the data-root env var.
- `notebooks/` run in numeric order; `exploratory/` is the alignment-quality side work.
- `results/tables/` (csv, tracked) and `results/figures/` (png) are the only outputs in git.
- `tests/` runs on synthetic fixtures built in `conftest.py`; no real data needed. `pytest` takes
  about three minutes because it fits models.
- `archive/` holds the previous generations (`scripts_pre_restructure/` is the script version of this
  pipeline, kept until the rerun has been compared; `Chromosome22_Oran_114/`, `whole_genome_Orangutan/`).
  `unmantained_code/` is explicitly unmaintained. Don't build on any of them or "clean them up"
  unless asked. Same-named files there differ from the package; diff before copying anything.

## Environment

One micromamba env, `Udenv` (`~/micromamba/envs/Udenv`, Python 3.14), runs everything including
downloads. Install: `micromamba create -f environment.yml`, then `pip install --no-deps -r
requirements-udenv.txt` and `pip install --no-deps -e .`. Pins that matter and why:

- `cogent3==2026.10.22`, `scinexus==2026.9.12`: the newer scinexus breaks `eti alignments`
  (`ensembl_tui/cli.py:633`). cogent3 declares a newer scinexus; `pip check` complains; expected.
- `duckdb<2`: duckdb 2.0 removed `fetchone()`, which `eti homologs` uses.
- `cogent3-h5seqs==0.9.1` while `ensembl_tui==0.7.9` pins 0.7.3: reading the existing `install/`
  works; writing a new install with 0.9.1 is unverified (fallback in `docs/pipeline_downloaddata.md`).
- `numba==0.68.1rc2` is the Python 3.14 build; a plain `numba` requirement will not resolve it.

The composable-app machinery is `scinexus` (`from scinexus.composable import define_app`);
`cogent3.app.composable` no longer exists. Fitted-model pickles are version-sensitive: keep one env
across stages 4 and 5. `environment.yml` lists only python and pip on purpose; the pip pins live in
`requirements-udenv.txt` (exact) and `pyproject.toml` (ranges, with comments).

Do not download anything from Ensembl or any other data repository, not even small test sets.

## Running the pipeline

```
export DATA_HUMCHIMPORANGOR114=/path/to/HumChimpOranGor_10_114
ud align-cds      --chrm $chr --submodel singlent   # stage 1: cds only
ud filter-introns --chrm $chr --submodel singlent   # stage 2: introns5UTR introns3UTR introns_nonUTR
ud filter-noncds  --chrm $chr --submodel singlent   # stage 3: intergenicAR intronsAR distalIG proximal5IG proximal3IG
ud human-seq      --chrm $chr                       # all 9 regions -> human_seq.fa
ud fit            --chrm $chr --submodel singlent   # stage 4: fit + pickle, all 9 regions
ud constraint                                       # stage 5: both submodels, chr 1-22 + X
ud quality <masked-databank|masked-concat|filter-spscore|fit-spscore> --chrm $chr [--window 60 --threshold 0.5]
```

`-chrm`/`-submodel` single-dash forms are accepted aliases. Repeat stages 1-4 with
`--submodel trinuc`; the submodel only changes `omit_degenerates(motif_length=1|3)` and the output
label. There is no motif-length option on `human-seq`: it always truncates to whole trinucleotides.

| Stage | Writes under `$DATA_HUMCHIMPORANGOR114/<region>/chrm<chr>/` |
| --- | --- |
| 1-3 | `singlent_filtered.fa` + `_alnstat.txt`, or `trinucleotide_filtered.fa` + `_alnstat.txt` |
| human-seq | `human_seq.fa` + `human_seq_alnstat.txt` |
| 4 | `sm_output/singlent_lh.pickle` or `sm_output/trinucleotide_lh.pickle` (`trinucleotide`, not `trinuc`) |
| 5 | `results/tables/{singlent,trinuc}_ENS.csv` and `_constraint.csv` |
| quality | `nonmasked_alldata_chrm<chr>/`, `singlent_nomasked.fa`, `quality_aligned_chrm<chr>/trinucleotide_filtered_window<W>_qualitythreshold<T>.fa` (+ `sm_output/`) |

Raw `eti` output in `<region>/alldata_chrm<chr>/` can be tarred and deleted once stages 1-3 and
`human-seq` have run.

### Notebook chain

1. `01_location_inter_intragenic.ipynb` runs *before* the intergenic downloads: it turns
   `eti dump-genes` output into `$DATA/intergenic_coordinates/chrom<chr>_{intergenic,distal_IG_coordinates,proximal5_IG_coordinates,proximal3_IG_coordinates}.tsv`.
   It is the one notebook that writes into the data root, by design.
2. `ud constraint` -> `02_summarystats_constraint.ipynb` (writes `*_constraint_mean_sem.csv`).
3. `ud human-seq` -> `03_estimating_mutrates.ipynb` (Oman et al. 2022 trinucleotide mutability,
   `mutability_motiflength3.csv`) -> `04_summarystats_mutability.ipynb` (normalises to 1.34e-8 per
   site; `average_mutrates_motiflength3_mean_sem.csv`, `UT_perregion_motiflength3.csv`, `length_perregion.csv`).
4. `05_Ud_estimation.ipynb` combines the mean rows of both `*_mean_sem.csv` files with the lengths.

Notebooks use the `udenv` kernel, read and write `results/tables` through `config.results_dir`, and
are committed with outputs (nbstripout runs with `--keep-output`; it strips execution counts only).

### SLURM

Submit from the repo root: `scripts/submit.sh scripts/slurm/04_fit.sh`. Logs go to
`logs/<jobname>_<jobid>_<task>.{out,err}` (gitignored). The download wrappers print a mask report
after each `eti alignments`; a chromosome where almost every fragment lacks `?` was downloaded
without masking and must be redone. `AttributeError: 'NotCompleted' object has no attribute
'write'` in a `.err` means an input alignment had length 0. `FileNotFoundError` for `cds/chrmY` in
stage 4 is expected; stage 5 skips Y.

## Things that bite

- **cds files split sequence names on `-`, every other region splits on `:`.** Using the wrong renamer
  silently drops every sequence, because each renamer ends with `take_seqs` on the renamed set.
- `stages/align_cds.py` pins `numpy.seterr(divide="ignore")` on purpose: cogent3 submodules disagree
  on the numpy error state at import time and `progressive_align` legitimately takes `log(0)`.
- `estimating_ud/__init__.py` imports `models` so the two trinuc models are registered; registration
  is guarded, so a reload does not raise.
- Stage 5 pads the three stop codons into the trinuc cds motif probabilities (1e-12 pseudocount)
  because `GNC_CpG_ss` has no stop states but the intergenicAR matrix does.
- The `scales["CpG"]` entry in both trinuc models is bound to `omega`, not to the CpG predicate
  (see `docs/todo.md`). Stage 5 never reads the scales, so leave it unless asked.
- The UTR samplers cut at the first/last `?` in the Human sequence, where both cds and repeats are
  masked, so "UTR" regions are a few hundred bp at each transcript end (`docs/todo.md` 0.1).
- `.gitignore` ignores `*.csv` except under `results/tables/` and the archived `output_data/` folders.
