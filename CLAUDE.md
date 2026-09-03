# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

A research pipeline (not a distributable package) for estimating Ud, the deleterious mutation rate, from
primate whole-genome alignments (Human/Chimpanzee/Orangutan, or Human/Chimpanzee/Gorilla in
`Chromosome22_Gorilla`), built on `cogent3` and `ensembl_tui` (CLI: `eti`). There is no build system, test
suite, or linter — work happens through numbered Python scripts (run per chromosome/region, usually via SLURM)
and Jupyter notebooks, operating on multi-terabyte alignment data stored outside the repo.

## Repository layout — folders are not shared modules

The repo root holds only documentation (`pipeline_downloaddata.md`, `pipeline_filtering_modelfitting.md`,
`README.md`), `environment.yml`, and an `__init__.py`. **All executable code lives in one of four
self-contained pipeline folders**, each with its own `paths.py` and `libs.py`:

- `whole_genome_Orangutan/` — the current, canonical Human/Chimpanzee/Orangutan (Ensembl 115) whole-genome
  pipeline. The two root `pipeline_*.md` files document *this* folder's scripts (they were written before the
  code was moved here in commit `0a81ff4`, so they still say "run `python3 codon_aligner_1.py`" without
  naming the folder).
- `Chromosome22_Oran_114/` — chromosome-22 pilot on Ensembl 114 with its own copies of `pipeline_*.md`. Its
  scripts diverge substantially from `whole_genome_Orangutan/` (see below) and are mid-refactor.
- `Chromosome22_Oran_115/` — notebooks only (`constraint_singlent.ipynb`, `constraint_trinuc.ipynb`) plus
  `libs.py`/`paths.py`. No scripts.
- `Chromosome22_Gorilla/` — notebook-driven Gorilla pipeline with its own `pipeline.md`, older package pins,
  and a `libs.py` that renames `gorilla_gorilla` → `Gorilla` and adds two extra
  `human_seq_length_nonfiltered_*` apps.
- `unmantained_code/` is explicitly unmaintained (per commit `1cca652`) — don't build on it or "clean it up"
  unless asked.

**Never assume a same-named script or module behaves identically across folders — diff before editing
shared-looking logic.** `libs.py` is byte-identical between `whole_genome_Orangutan/`,
`Chromosome22_Oran_114/`, and `Chromosome22_Oran_115/`; everything else differs.

## Environment setup

`environment.yml` is the frozen snapshot of the working env and is named `UdChimpHumOran` — that is also the
env every `whole_genome_Orangutan/bash_*.sh` wrapper activates. The `pipeline_*.md` docs name *different*
envs, one per pipeline generation, and pin different package versions:

| Pipeline folder | Env name in its docs | Pins |
| --- | --- | --- |
| root docs / `whole_genome_Orangutan` | `Ensembl0.7.6` | `ensembl_tui==0.7.6`, `cogent3==2026.1.12a1`, `cogent3_h5seqs==0.7.3` |
| `Chromosome22_Oran_114` | `Ensembl0.7.7` | `ensembl_tui==0.7.7`, `cogent3==2026.7.30a0`, `cogent3_h5seqs==0.7.3` |
| `Chromosome22_Gorilla` | `Udestimation_env` | `cogent3==2025.7.10a5`, `cogent3_h5seqs==0.5.0`, `ensembl_tui==0.4.3` |

Always read the `pipeline*.md` in the folder you are working in for its own pins. Plus `pandas matplotlib`
(and `phylim`, which `libs.py` imports unconditionally).

**Before running anything**, edit the `paths.py` *of the folder you are running from* so its data-root
variable points at your alignment download location. The variable name differs per folder and scripts
reference it by name:

- `whole_genome_Orangutan/paths.py`, `Chromosome22_Oran_115/paths.py` → `DATA_HUMCHIMPORANG115`
- `Chromosome22_Oran_114/paths.py` → `DATA_HUMCHIMPORANGOR114`
- `Chromosome22_Gorilla/paths.py` → `DATA_HUMCHIMPGOR115`

## Running the pipeline

Read both root docs before changing filtering or model-fitting logic:
`pipeline_downloaddata.md` (installing alignments per genomic region via `eti`) and
`pipeline_filtering_modelfitting.md` (how the scripts chain together, plus the storage-reclamation recipe).

From inside `whole_genome_Orangutan/`, stages run per chromosome (`-chrm`, values `1`–`22`, `X`, `Y`) and/or
region (`-reg`), each writing back under `paths.DATA_HUMCHIMPORANG115/<region>/chrm<chr>/`. Every stage has a
single-nucleotide and a trinucleotide variant; the trinucleotide ones differ only in
`omit_degenerates(motif_length=3)` instead of `1` and in writing `trinucleotide_filtered.fa` /
`trinucleotide_filtered_alnstat.txt` instead of `filtered.fa` / `filtered_alnstat.txt`.

1. `codon_aligner_1.py -chrm $chr` (cds only) — codon-aligns, trims stop codons → `cds/chrm$chr/filtered.fa`.
   Trinuc counterpart: `codon_aligner_trinucs_1.py`.
2. `filtering_gaps_introns_2.py -chrm $chr` — splits intron alignments into `introns5UTR`/`introns3UTR`/
   `introns_nonUTR` (three hardcoded blocks in one run) and strips gapped/degenerate sites.
   Trinuc counterpart: `filtering_trinucs_gaps_introns_2.py`.
3. `filtering_gaps_nointrons_3.py -reg <region> -chrm $chr` — same stripping for `intergenicAR`, `intronsAR`,
   `distalIG`, `proximal5IG`, `proximal3IG`. Trinuc counterpart: `filtering_trinucs_gaps_nointrons_3.py`.
4. `fittinglh_singlentmodel_4.py -reg <region> -chrm $chr` — fits a `GN` model with
   `time_het="max"`, `discrete_edges=["Orangutan"]`, reading `filtered.fa` → `sm_output/singlent_lh.pickle`
   and `sm_output/singlent_alnstat.txt`. `fittinglh_trinucmodel_4.py` reads `trinucleotide_filtered.fa` and
   fits the trinucleotide CpG models from `trinuc_models.py` → `sm_output/trinuc_lh.pickle`,
   `sm_output/trinuc_alnstat.txt`.
5. `meassuring_trinuc_constraint_5.py` — takes no arguments; sweeps every region × chromosome (chr 1–22 and X,
   **not** Y), compares each region's ENS against the `intergenicAR` neutral rate, and writes
   `output_data/trinuc_ENS.csv`, `output_data/trinuc_constraint.csv`,
   `output_data/trinuc_constraint_mean_sem.csv`. Those paths are **relative to the CWD and the script does not
   create the directory** — `mkdir -p output_data` first. `*.csv` is gitignored, so these outputs are never
   committed.
6. `Ud_estimation.ipynb` / `Ud_estimation_approx.ipynb` — final Ud calculation; each requires the constraint
   output above and `numb_sites_region.ipynb` to have already been run.

### SLURM

Each stage has a matching wrapper: `bash_codonaliger1.sh`, `bash_intron_filtering2.sh`, `bash_filtering3.sh`,
`bash_singlent_submodel4.sh`, and the `bash_trinuc_*` variants. They are array jobs that decompose
`SLURM_ARRAY_TASK_ID` into a region × chromosome pair (`1-24` chromosome-only, `1-120` for 5 regions × 24
chromosomes, `1-216` for 9 regions × 24), activate `UdChimpHumOran`, `mkdir -p logs_*/`, and invoke the script
from the CWD — so **submit them from inside `whole_genome_Orangutan/`**.

When a job fails, check the `.err` files under `logs_*/`. `AttributeError: 'NotCompleted' object has no
attribute 'write'` almost always just means the input alignment had length 0 and can be ignored — use the
`find`/`grep` one-liner in `pipeline_filtering_modelfitting.md` to isolate real failures.

## Architecture

- `libs.py` — shared `cogent3` composable apps. Four species renamers: **cds files split sequence names on
  `-`, non-cds files split on `:`** — using the wrong one silently drops every sequence, since each renamer
  ends with `take_seqs(list(name_map.values()))`. Also the UTR5/UTR3/non-UTR intron samplers (which locate
  UTR boundaries by the first/last `?` mask character in the Human gapped sequence), site/motif counters, a
  `phylim` split-codon identifiability workaround, and the dinucleotide CpG model builders `GDN_CpG` /
  `GDN_CpG_ss`.
- `trinuc_models.py` — registers `GNC_CpG_ss` (codon states, for cds) and `GT_CpG_ss` (all trinucleotide
  states, for non-cds) as `cogent3` "codon" models, both with GN predicates, `omega`, strand-symmetric CpG
  deamination, and `scales={"CpG", "notCpG"}` — the two scales are what stage 5 reads via
  `get_scaled_lengths`. `modified_lf` rebuilds a likelihood function from another region's fitted parameter
  rules, which is how a region's ENS is recomputed under the `intergenicAR` rate matrix.
- Filtering scripts share one shape: `load_*` → rename species → region-specific transform →
  `omit_degenerates` → `concat` across per-gene fragments → write `filtered.fa` + `filtered_alnstat.txt`.
  Model-fitting scripts share another: load a `filtered.fa`, fit, pickle the `model_result`, write an alnstat
  file. Match these shapes rather than inventing a new one for a new region/script.
- Regions fall into two families with different loaders: `cds` starts as unaligned FASTA (from `eti homologs`)
  and is codon-aligned by the pipeline itself (`load_unaligned` + `progressive_align`); everything else
  (`introns*`, `intergenicAR`, `intronsAR`, `distalIG`, `proximal5IG`, `proximal3IG`) arrives pre-aligned
  from `eti alignments` and only needs gap/degenerate filtering (`load_aligned`).
- Data layout on disk: `$region/alldata_chrm$chr` holds raw `eti` output (one file per contig, terabytes);
  `$region/chrm$chr` holds the pipeline's concatenated filtered alignment and `sm_output/`. The `alldata_*`
  folders can be tarred and deleted once filtering has run.
- Species are renamed early in every pipeline to `Human`/`Chimpanzee`/`Orangutan` (or `Gorilla`) — don't rely
  on raw Ensembl species IDs (`homo_sapiens`, `pan_troglodytes`, `pongo_abelii`) downstream of the renamer.

## Known rough edges (verify before trusting a script)

The move of the scripts into `whole_genome_Orangutan/` left the import style inconsistent, and
`Chromosome22_Oran_114` is mid-refactor. Check these before assuming a script runs:

- **Mixed import style within `whole_genome_Orangutan/`.** Most scripts do `import paths` / `import libs`
  (must run with that folder as CWD), but `codon_aligner_trinucs_1.py` does
  `import whole_genome_Orangutan.paths as paths` (must run from the repo root). Its wrapper
  `bash_trinuc_codonaliger1.sh` invokes it from the folder, so the two disagree.
- **`Chromosome22_Oran_114` scripts mostly read the *other* folder's data root.** All of them except
  `codon_aligner_1.py` import `whole_genome_Orangutan.paths` and use `DATA_HUMCHIMPORANG115`, so the local
  `DATA_HUMCHIMPORANGOR114` in `Chromosome22_Oran_114/paths.py` is only honoured by `codon_aligner_1.py`.
- **Different CLI contracts in `Chromosome22_Oran_114`.** Its scripts take `-submodel {singlent,trinuc}`
  (selecting `omit_degenerates` motif length 1 vs 3 at runtime) and loop over all regions internally, instead
  of the root pipeline's fixed-model, one-region-per-invocation `-reg`/`-chrm` design. Several also hardcode
  `chrm22`. `fittinglh_singlentmodel_4.py` still reads `args.chromosome` although it no longer defines
  `-chrm`, and `filtering_gaps_nointrons_3.py` tests `args.substitutionmodel == False` against a string, so it
  always takes the trinucleotide branch.
- `Chromosome22_Oran_114/trinuc_models.py` appends an experimental block after `modified_lf`
  (`modified_lf_cds`, `modified_lf_Ulisesmodified`, and Q/f0/length parameter helpers) that the root version
  does not have; its own comments mark it as possibly disposable.
