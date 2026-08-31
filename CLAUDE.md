# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

A research pipeline (not a distributable package) for estimating Ud, the deleterious mutation rate, from
primate whole-genome alignments (Human/Chimpanzee/Orangutan, or Human/Chimpanzee/Gorilla in the `_Gorilla`
folder), built on `cogent3` and `ensembl_tui` (CLI: `eti`). There is no build system, test suite, or linter —
work happens through numbered Python scripts (meant to run per chromosome/region, often via SLURM) and Jupyter
notebooks, operating on multi-terabyte alignment data stored outside the repo.

## Environment setup

```
conda create -n UdChimpHumOran python=3.13
conda activate UdChimpHumOran
pip install "ensembl_tui==0.7.6"
pip install "cogent3==2026.1.12a1" "cogent3_h5seqs==0.7.3" -U
pip install pandas matplotlib phylim
```

- Package versions are pinned per pipeline; check the `pipeline_downloaddata.md`/`pipeline.md` in the folder
  you're working in — the `Chromosome22_Oran_114` pilot uses different pins (`cogent3==2026.7.30a0`,
  `ensembl_tui==0.7.7`) and a differently named env (`UdEst`).
- Full frozen environment snapshot: `environment.yml`.
- **Before running anything**, edit `paths.py` (root, or the per-folder copy) so its data-root variable points
  at your own alignment download location — every script imports it directly.
- There is no test runner, linter, or build step in this repo.

## Running the pipeline

The root pipeline is documented across two files — read them before changing filtering/model-fitting logic:
- `pipeline_downloaddata.md` — installing alignments per genomic region via `eti`
- `pipeline_filtering_modelfitting.md` — how the filtering/model-fitting scripts below chain together

Numbered scripts run in stage order per chromosome (`-chrm`) and/or region (`-reg`), each writing outputs back
under `paths.DATA_HUMCHIMPORANG115/<region>/chrm<chr>/...`:

1. `codon_aligner_1.py -chrm $chr` (cds only) — codon-aligns, trims stop codons → `cds/chrm$chr/filtered.fa`.
   `codon_aligner_trinucs_1.py` is the trinucleotide-model counterpart → `trinucleotide_filtered.fa`.
2. `filtering_gaps_introns_2.py -chrm $chr` — splits intron alignments into `introns5UTR`/`introns3UTR`/
   `introns_nonUTR`, strips gapped/degenerate sites → `filtered.fa` per sub-region.
   `filtering_trinucs_gaps_introns_2.py` is the trinuc counterpart.
3. `filtering_gaps_nointrons_3.py -reg <region> -chrm $chr` — same gap/degenerate stripping for
   `intergenicAR`, `intronsAR`, `distalIG`, `proximal5IG`, `proximal3IG`.
   `filtering_trinucs_gaps_nointrons_3.py` is the trinuc counterpart.
4. `fittinglh_singlentmodel_4.py -reg <region> -chrm $chr` — fits a single-nucleotide GN substitution model,
   pickles the result to `sm_output/singlent_lh.pickle`. `fittinglh_trinucmodel_4.py` fits the trinucleotide
   CpG model instead (`sm_output/trinuc_lh.pickle`), using models from `trinuc_models.py`.
5. `meassuring_trinuc_constraint_5.py` — sweeps every region/chromosome, compares each region's ENS against
   the `intergenicAR` neutral rate to compute constraint, writes `output_data/trinuc_ENS.csv`,
   `output_data/trinuc_constraint.csv`, `output_data/trinuc_constraint_mean_sem.csv`.
6. `Ud_estimation.ipynb` / `Ud_estimation_approx.ipynb` — final Ud calculation; each requires the constraint
   output above and `numb_sites_region.ipynb` to have already been run.

On a SLURM cluster, use the matching `bash_*.sh` wrapper (`bash_codonaliger1.sh`, `bash_intron_filtering2.sh`,
`bash_filtering3.sh`, `bash_singlent_submodel4.sh`, and the `bash_trinuc_*` variants) — these submit array jobs
sweeping all region×chromosome combinations. When a job fails, check the `.err` files under `logs*/`;
`AttributeError: 'NotCompleted' object has no attribute 'write'` almost always just means the input alignment
had length 0 and can be ignored — use the `find`/`grep` one-liner in `pipeline_filtering_modelfitting.md` to
isolate real failures.

## Architecture

- `paths.py` — single source of truth for the external data root (e.g. `DATA_HUMCHIMPORANG115`); every script
  imports it. This is the one thing a new user must edit (per `README.md`).
- `libs.py` / `__init__.py` — shared `cogent3` composable apps: species renamers (cds files split names on
  `-`, non-cds files split on `:` — use the matching renamer), UTR5/UTR3/non-UTR intron samplers, and
  nucleotide/dinucleotide substitution model builders (`GDN_CpG`, `GDN_CpG_ss`).
- `trinuc_models.py` — trinucleotide/codon substitution model registrations (`GNC_CpG_ss` for cds, `GT_CpG_ss`
  for non-cds) plus `modified_lf`, used to recompute a region's ENS under another region's fitted rate matrix.
- Filtering scripts share one shape: `load_*` → rename species → region-specific transform →
  `omit_degenerates` → `concat` across per-gene fragments → write `filtered.fa` + `filtered_alnstat.txt`.
  Model-fitting scripts share another: load a `filtered.fa`, fit, pickle the `model_result`, write an alnstat
  file. Match these shapes rather than inventing a new one for a new region/script.
- Regions fall into two families with different loaders: `cds` starts as unaligned FASTA and is codon-aligned
  by the pipeline itself (`load_unaligned` + `progressive_align`); everything else (`introns*`, `intergenicAR`,
  `intronsAR`, `distalIG`, `proximal5IG`, `proximal3IG`) arrives pre-aligned from `eti alignments` and only
  needs gap/degenerate filtering (`load_aligned`).
- Species are renamed early in every pipeline to `Human`/`Chimpanzee`/`Orangutan` (or `Gorilla`) — don't rely
  on raw Ensembl species IDs (`homo_sapiens`, `pongo_abelii`, ...) downstream of the renamer step.

## Repository layout — folders are not shared modules

- Root directory holds the current Human/Chimpanzee/Orangutan (Ensembl release 115) pipeline described above.
- `Chromosome22_Oran_114/`, `Chromosome22_Oran_115/`, `Chromosome22_Gorilla/` are self-contained prior/parallel
  pipeline runs, each with its own `paths.py`, `libs.py`, `pipeline.md`, and often its own copy of scripts with
  real behavioral divergences from the root version — e.g. `Chromosome22_Oran_114/filtering_gaps_introns_2.py`
  takes a `-submodel {singlent,trinuc}` flag and loops all three intron sub-regions in one run, instead of the
  root's fixed-single-nucleotide, hardcoded-triple-block version. **Never assume a same-named script behaves
  identically across folders — diff before editing shared-looking logic.**
- `unmantained_code/` is explicitly unmaintained (per commit history) — don't build on it or "clean it up"
  unless asked.
- `output_data/` holds the committed CSV outputs of `meassuring_trinuc_constraint_5.py`. Per-folder
  `Output_data/`/`output_data/` directories hold the notebook-produced equivalents for that folder's pipeline.
