# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

A research pipeline (not a distributable package) for estimating Ud, the deleterious mutation rate, from
Human/Chimpanzee/Orangutan whole-genome alignments (Ensembl 114, downloaded with `ensembl_tui`, CLI `eti`)
using `cogent3` substitution models. There is no build system, test suite, or linter. Work happens through
numbered Python scripts run per chromosome (usually as SLURM array jobs) and Jupyter notebooks, operating on
multi-terabyte alignment data stored outside the repo.

## Repository layout

**The repo root is the canonical, current pipeline.** Everything else is an older generation:

- Root — scripts `codon_aligner_1.py` … `estimating_constraint_5.py`, `get_human_sequence.py`, `libs.py`,
  `trinuc_models.py`, `paths.py`, the `bash_*.sh` SLURM wrappers, the analysis notebooks, and the three docs
  `pipeline_downloaddata.md`, `pipeline_filtering.md`, `pipeline_constraint_estimates.md`. Read those docs
  before changing download, filtering, or model-fitting logic.
- `Chromosome22_Oran_114/` — the chromosome-22 pilot the root scripts were generalised from. Same script
  names, but they hardcode `chrm22`, take no `-chrm`, its `paths.py` points at a laptop path, its docs name env
  `Ensembl0.7.7`, and its `trinuc_models.py` carries an experimental block after `modified_lf`
  (`modified_lf_cds`, `modified_lf_Ulisesmodified`, Q/f0/length helpers) that its own comments mark as
  disposable. `Chromosome22_Oran_114/test/` holds scratch notebooks.
- `whole_genome_Orangutan/` — the previous generation on Ensembl 115 (`DATA_HUMCHIMPORANG115`, env
  `UdChimpHumOran`, one-region-per-invocation `-reg`/`-chrm` CLI, separate `*_trinucs_*` scripts). Its
  `libs.py` imports `cogent3.app.composable`, which no longer exists in the current cogent3, so it will not
  import under `Ensembl0.7.9`. Superseded; don't build on it.
- `unmantained_code/` — explicitly unmaintained (commit `1cca652`), including the older Gorilla pipelines.
  Don't build on it or "clean it up" unless asked.

**Never assume a same-named script behaves identically across folders — diff before editing.** `libs.py`
differs between all three script folders (root has the `gethumanseq_*` apps and a `motif_length`-aware
`number_of_motifs`).

## Environment

Use the conda env `Ensembl0.7.9` (at `~/.conda/envs/Ensembl0.7.9`; it is what every root `bash_*.sh`
activates and what the notebooks' kernel is named). Pins that matter: `ensembl_tui==0.7.9`,
`cogent3==2026.9.10`, `cogent3-h5seqs==0.7.3`, plus `scinexus`, `pandas`, `matplotlib`, `seaborn`.

- In this cogent3, the composable-app machinery lives in the `scinexus` package:
  `from scinexus.composable import define_app` (root `libs.py`) is correct, and `cogent3.app.composable` is
  gone. `cogent3.app.typing` still exists.
- `environment.yml` at the root is a **stale** snapshot of the older `UdChimpHumOran` env (cogent3 2026.1,
  ensembl_tui 0.7.6, phylim) — don't recreate the env from it and don't trust its pins.
- Fitted-model pickles are sensitive to library versions; keep one env across stages 4 and 5.

Before running anything, set `DATA_HUMCHIMPORANGOR114` in the root `paths.py` to the alignment download
location (the README still calls this variable `DATA_APES114`; that name is stale). Every root script and
notebook reads paths through that variable.

## Running the pipeline

All commands run from the repo root with `Ensembl0.7.9` active. Scripts take `-chrm` (one of `1`–`22`, `X`,
`Y`) and `-submodel {singlent,trinuc}`; the submodel only changes `omit_degenerates(motif_length=1|3)` and the
output filename. Each script loops over its own hardcoded region list internally — there is no `-reg` flag.

```
python3 codon_aligner_1.py          -chrm $chr -submodel singlent   # cds only: codon-align, trim stops
python3 filtering_gaps_introns_2.py -chrm $chr -submodel singlent   # introns5UTR / introns3UTR / introns_nonUTR
python3 filtering_gaps_nointrons_3.py -chrm $chr -submodel singlent # intergenicAR intronsAR distalIG proximal5IG proximal3IG
python3 get_human_sequence.py       -chrm $chr -mutmotif 3          # all 9 regions: concatenated Human sequence
python3 fittinglh_sm_4.py           -chrm $chr -submodel singlent   # all 9 regions: fit + pickle
mkdir -p output_data && python3 estimating_constraint_5.py          # no args; both submodels, chr 1-22 + X
```

Repeat stages 1–4 with `-submodel trinuc`. Outputs land under `$DATA_HUMCHIMPORANGOR114/<region>/chrm<chr>/`:

| Stage | Writes |
| --- | --- |
| 1–3 | `singlent_filtered.fa` + `singlent_filtered_alnstat.txt`, or `trinucleotide_filtered.fa` + `trinucleotide_filtered_alnstat.txt` |
| `get_human_sequence.py` | `human_seq_motiflength<N>.fa` + `_alnstat.txt` (N = `-mutmotif`) |
| 4 | `sm_output/singlent_lh.pickle` or `sm_output/trinucleotide_lh.pickle` (note: `trinucleotide`, not `trinuc`) |
| 5 | `output_data/{singlent,trinuc}_ENS.csv`, `_constraint.csv`, `_constraint_mean_sem.csv` (CWD-relative; `*.csv` is gitignored) |

Raw `eti` output lives in `<region>/alldata_chrm<chr>/` (one file per contig, terabytes) and can be tarred
and deleted once stages 1–3 and `get_human_sequence.py` have run — see `pipeline_filtering.md`.

### Notebook chain after the scripts

- `location_inter_intragenic.ipynb` runs *before* the intergenic downloads: it turns `eti dump-genes` output
  into `intergenic_coordinates/chrom<chr>_{distalIG,proximal5IG,proximal3IG}_coordinates.tsv` under the data root.
- Constraint: `estimating_constraint_5.py` → `plot_constraint.ipynb` (plots, and rewrites the
  `*_constraint_mean_sem.csv` files).
- Mutation rate: `get_human_sequence.py -mutmotif 3` → `estimating_mutrates.ipynb` (applies the Oman et al.
  2022 trinucleotide mutability model, embedded as a dict, writes `output_data/mutability_motiflength3.csv`)
  → `plot_mutability.ipynb` (normalises to a per-site rate of 1.34e-8, writes
  `output_data/average_mutrates_motiflength3_mean_sem.csv` and `output_data/length_perregion.csv`).
- `Ud_estimation.ipynb` reads `trinuc_constraint_mean_sem.csv`, `length_perregion.csv`, and
  `average_mutrates_motiflength3_mean_sem.csv` and computes Ud. It needs both chains above.
- The notebooks hardcode `motif_length = 3`, and the human-sequence filenames carry the motif length, so
  `get_human_sequence.py` must be run with `-mutmotif 3` for them to find their input.
- `mutation_rate_perregion.ipynb` is partly stale (its later cells reference `paths.DATA_APES114`, which no
  longer exists); `estimating_mutrates.ipynb` + `plot_mutability.ipynb` supersede it.

### SLURM

`bash_codonaliger1.sh`, `bash_intron_filtering2.sh`, `bash_filtering3.sh`, `bash_get_humanseq.sh`,
`bash_submodel4.sh` are `--array=1-24` jobs (task ID → chromosome), each running the `singlent` then `trinuc`
invocation of one script. They `mkdir -p logs_*/`, activate `Ensembl0.7.9`, and call the script from the CWD,
so **submit them from the repo root** (`sbatch bash_submodel4.sh`). `logs*` is gitignored.

When a task fails, read the `.err` under `logs_*/`. `AttributeError: 'NotCompleted' object has no attribute
'write'` means an input alignment had length 0 and can be ignored; use the `find`/`grep -FL` one-liner in
`pipeline_filtering.md` to list the other failures. `FileNotFoundError` for `cds/chrmY/*_filtered.fa` in
stage 4 is expected (no cds alignment for Y), and stage 5 skips Y entirely.

## Architecture

- `libs.py` — `cogent3` composable apps. Four species renamers map `homo_sapiens`/`pan_troglodytes`/
  `pongo_abelii` → `Human`/`Chimpanzee`/`Orangutan`: **cds files split sequence names on `-`, non-cds files
  split on `:`**. Using the wrong one silently drops every sequence, since each ends with
  `take_seqs(list(name_map.values()))`. Also `gethumanseq_cds_unaligned` / `gethumanseq_noncds_aligned`
  (extract the Human sequence, truncated to whole motifs), the UTR5/UTR3/non-UTR intron samplers (which find
  UTR boundaries by the first/last `?` mask character in the Human gapped sequence), site/motif counters, and
  the dinucleotide CpG models `GDN_CpG` / `GDN_CpG_ss`.
- `trinuc_models.py` — registers `GNC_CpG_ss` (codon states, for cds) and `GT_CpG_ss` (all 64 trinucleotide
  states, for non-cds) via `@register_model("codon")` as an *import side effect*: `import trinuc_models`
  before `get_app("model", "GT_CpG_ss", ...)` or the name is unknown. Both have GN predicates, `omega`,
  strand-symmetric CpG deamination, and `scales={"CpG", "notCpG"}`. `modified_lf` rebuilds a likelihood
  function from another fit's parameter rules.
- Model choices in `fittinglh_sm_4.py`: single-nucleotide is `GN` with `time_het="max"` and
  `discrete_edges=["Orangutan"]`; trinucleotide cds fits `GNC_CpG_ss` with all params independent, non-cds
  fits `GT_CpG_ss` with `omega` fixed at 1.
- `estimating_constraint_5.py` computes each region's ENS (expected number of substitutions on the Human
  edge) under its own fitted Q and under the `intergenicAR` Q of the same chromosome, then
  `constraint = (ENS_IGAR − ENS)/ENS_IGAR`. `intergenicAR` is therefore not in its `REGIONS` list; it appears
  in the ENS csv as `IGAR`. For trinuc cds it pads the three stop codons into the motif probabilities with a
  1e-12 pseudocount because `GNC_CpG_ss` has no stop-codon states but the IGAR matrix does.
- Two region families with different loaders: `cds` arrives as unaligned FASTA from `eti homologs` and is
  codon-aligned by the pipeline (`load_unaligned` + `progressive_align`, guide tree
  `(Human:0.06,Chimpanzee:0.06,Orangutan:0.14)`); everything else arrives pre-aligned from `eti alignments`
  (`load_aligned`) and only needs gap/degenerate filtering. The three intron regions all read from
  `introns/alldata_chrm<chr>` and differ only in the sampler app.
- Every filtering script has one shape: `load_*` → renamer → region transform → `omit_degenerates` →
  `concat` across per-gene fragments → write `.fa` + `_alnstat.txt`. Match this shape for new regions.
- `codon_aligner_1.py` pins `numpy.seterr(divide="ignore")` on purpose: cogent3 submodules disagree on the
  numpy error state at import time, and `progressive_align` legitimately takes `log(0)`. Don't remove it.
