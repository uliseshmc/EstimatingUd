# Inconsistencies found in the root-folder pipeline (2026-09-24)

Read-only audit of the scripts, SLURM wrappers, notebooks, and `output_data/` CSVs in the repo root.
Nothing was changed. Items are grouped by how much they affect the numbers currently on disk.

## 0. HIGHEST PRIORITY: repeat masking failed at download time for 5 of 6 non-cds regions (found 2026-09-24)

Checked the raw `eti alignments` output under `$DATA/<region>/alldata_chrm*/` by counting `?` in the
homo_sapiens sequence of a random sample of files per chromosome.

| Region | Download flag | Unmasked share of human bases, chr 1/5/12/19/21/X | chr22 |
| --- | --- | --- | --- |
| introns | `--mask cds_allAR --coord_names` | 40-54% (masking worked) | 49% |
| intronsAR | `--mask_shadow ancestralrepeats --coord_names` | **100% (no `?` in any file)** | 45% |
| intergenicAR | `--mask_shadow ancestralrepeats --ref_coords` | **100%** | 55% |
| distalIG | `--mask allAR --ref_coords` | **100%** | 33% |
| proximal5IG | `--mask allAR --ref_coords` | **100%** | 40% |
| proximal3IG | `--mask allAR --ref_coords` | **100%** | 43% |

Only chromosome 22 has working masks for every region (it was the test chromosome; see the "Bug caution"
comments in `pipeline_downloaddata.md`, which describe exactly this symptom: "files do not have any masked
positions (No question marks)"). Every `--mask_shadow` run and every `--mask ... --ref_coords` run on the
other chromosomes produced unmasked sequence.

Consequences for the numbers currently on disk:
- `intronsAR` is the entire aligned protein-coding gene body (CDS + UTR exons + all introns + all repeats),
  not intronic ancestral repeats. Its total, 1341 Mb, matches the merged protein-coding gene span computed
  from `homo_sapiens-114-gene_metadata.tsv` (1306 Mb, 45% of the non-N genome). Real intronic TE content
  should be roughly half of that.
- `intergenicAR` (745 Mb) is all aligned intergenic sequence, and `distalIG + proximal5IG + proximal3IG`
  (756 Mb) is the same sequence again without repeat removal. That is why distalIG constraint is ~0 (-0.02):
  it is being compared with itself.
- `cds` and `introns_*` are subsets of `intronsAR`, and the intergenic regions are subsets of the neutral
  reference, so `Ud_estimation.ipynb` double counts and the regions sum to 121% of the genome.
- The `intronsAR` constraint (0.13) is the constraint of whole gene bodies (incl. CDS) vs. all intergenic DNA.

Fix: re-run the five affected `eti alignments` downloads (`bashcommands/bash_sampleintronsAR.sh`,
`bash_sampleIGAR.sh`, `bash_sampledistalIG.sh`, `bash_sample5IG.sh`, `bash_sample3IG.sh`) once the
`ensembl_tui` masking bug with `--mask_shadow` / `--ref_coords` is understood (chr22 shows it can work; find
out what differed for that run: eti version, coordinate file format, or the `--mask_ref` flag used in the
test commands). Verify with `grep -c '?'` on the human sequences before running any filtering. Then re-run
stages 2-5 and the mutability notebook for those regions.

### 0.1 The UTR split does not delimit UTRs
In the `introns` alignments both CDS and repeats are masked with the same `?` character, in hundreds of
short runs per gene. `sample_UTR5`/`sample_UTR3`/`removeUTRs_fromintrons` (`libs.py:126-168`) cut at the
first/last `?`, i.e. at the nearest masked repeat or CDS from the transcript ends, typically a few hundred bp
in. Hence introns5UTR and introns3UTR are both ~6.7 Mb genome-wide (~300 bp per gene) and nearly equal, while
real 5'UTR exons total ~6-9 Mb and 3'UTR exons ~28-35 Mb (Ensembl canonical transcripts). To sample UTRs
properly the CDS mask must be distinguishable from the repeat mask (e.g. two downloads, or use exon
coordinates from `eti dump-genes`).

## 1. Affects numbers currently in `output_data/`

### 1.1 `estimating_constraint_5.py` accumulates rows across submodels
`row_data_ENS` and `row_data_constraint` are created once at `estimating_constraint_5.py:104-105`, outside
the `for submodel in SUBMODELS` loop, and never reset. `SUBMODELS = ["trinuc", "singlent"]`, so the trinuc
pass writes correct files and the singlent pass then appends onto the same lists.

Verified: the first 184 rows of `singlent_constraint.csv` are exactly `trinuc_constraint.csv`, and the first
391 rows of `singlent_ENS.csv` are exactly `trinuc_ENS.csv`. `singlent_constraint_mean_sem.csv` therefore
averages trinuc and singlent constraints together.

| File | Rows per region | Expected |
| --- | --- | --- |
| singlent_constraint.csv | 46 | 23 |
| singlent_ENS.csv | 92 (IGAR 46) | 46 (IGAR 23) |
| trinuc_constraint.csv | 23 | 23 (OK) |
| trinuc_ENS.csv | 46 | 46 (OK, selfQ 1 and 0 by design) |

Fix: move the two list initialisations inside the submodel loop.

### 1.2 `mutability_motiflength3.csv` has introns_nonUTR duplicated
Every chromosome of introns_nonUTR appears twice with identical values (46 rows instead of 23; 207 rows total
instead of 184). Region order in the file is cds, introns_nonUTR, introns3UTR, introns5UTR, introns_nonUTR,
intronsAR, ..., so `REGIONS` in `estimating_mutrates.ipynb` listed introns_nonUTR twice when the cell was
last run. The current notebook code lists it once and would produce a clean file on rerun.

Downstream: `length_perregion.csv` has a doubled introns_nonUTR site count (and deflated fractions for every
other region), `average_mutrates_motiflength3_mean_sem.csv` understates its SEM, the normalising constant is
biased toward intronic mutability, and `Ud_estimation.ipynb`'s Ud = 13.15 was computed from the doubled count.

### 1.3 Mean/SEM CSV format differs between script and notebooks
`estimating_constraint_5.py:173` writes `*_constraint_mean_sem.csv` with `index=False`, which drops the
Region labels (single column named `0`). `plot_constraint.ipynb` rewrites the same files with the index.
`Ud_estimation.ipynb` expects the notebook format (filters column 0 == "mean", takes columns 1:3).

Both mean_sem files currently on disk are in the script format (no Region column), so re-running
`Ud_estimation.ipynb` now would fail or read the wrong columns. Pick one format.

### 1.4 Two different definitions of "number of sites"
- `*_ENS.csv` `aln_length`: filtered alignment length (all three species ungapped and non-degenerate).
- `mutability_motiflength3.csv` `numb_sites` and `length_perregion.csv`: raw human sequence length from
  `get_human_sequence.py`, with no gap/degenerate filtering; for cds also no codon alignment or stop trimming.

Example, cds chr22: 661,098 (human seq) vs 536,267 (filtered alignment). `Ud_estimation.ipynb` multiplies a
constraint estimated on the first set of sites by a site count from the second. Possibly intentional, but
the two are not interchangeable and the choice should be documented.

### 1.5 Mutation-rate means are unweighted across chromosomes
`get_normalizer_mut_constant` in `plot_mutability.ipynb` weights each region x chromosome row by its
`numb_sites`, so the site-weighted mean of the per-row rate equals `persite_mu` exactly (verified: 1.3400e-8).
`stats_mutrate` then takes a plain `groupby("Region").mean()` across chromosomes, and cell 9 weights that
unweighted mean by region size, giving 1.3612e-8. Small high-mutability chromosomes (chr 19, 22, 16, 17) are
over-weighted. The plotted relative rates and the per-region rates fed into Ud have the same issue.

Fix: use a site-weighted mean across chromosomes for the per-region rate, or accept the discrepancy and say so.

## 2. Bugs that do not affect current numbers but will bite

### 2.1 The "CpG" scale in both trinucleotide models is bound to omega
`trinuc_models.py:63-64` and `:85-87` build `ssym_preds = gn_preds + cpg_preds + omega` and then set
`cpgdecay = ssym_preds[-1]`, which is omega, not the CpG predicate. Verified by building `GT_CpG_ss()` and
inspecting `scale_masks["CpG"]`: it selects 438 changes such as TTT->TTA, only 46 of which have a CG source.
So `scales["CpG"]` is really "nonsynonymous" and `scales["notCpG"]` is "synonymous".

`estimating_constraint_5.py` uses `expected_number_subs` and never reads the scales, so current constraint
values are unaffected. Anything calling `get_scaled_lengths` (e.g. stage 5 in `whole_genome_Orangutan/`) gets
mislabeled quantities. The identical line exists in `whole_genome_Orangutan/trinuc_models.py`.

Fix: `cpgdecay = ssym_preds[-2]`, or name the predicate explicitly before building the list.

### 2.2 Swapped alnstat file names
`estimating_constraint_5.py:41-44`: singlent reads `trinucleotide_filtered_alnstat.txt`, trinuc reads
`singlent_filtered_alnstat.txt`. Only the `aln_length` column of the ENS files is wrong.

### 2.3 `-mutmotif` argument to `get_human_sequence.py` is ignored
`get_humanseq_app` (`get_human_sequence.py:45-53`) hardcodes motif length 3 in every branch; the CLI value
only changes the output file name. The uncommitted edit to `bash_get_humanseq.sh:32` passes `-mutmotif 1`,
so a rerun would write `human_seq_motiflength1.fa` containing sequences truncated to multiples of 3, while
`estimating_mutrates.ipynb` keeps reading the motiflength3 file. `pipeline_filtering.md` says to use 3.

### 2.4 Intron alignments with no exon mask are handled three different ways
`libs.py:126-168`: when the Human gapped sequence has no `?`,
- `sample_UTR5` returns the whole alignment as 5' UTR (`aln[0:None]`),
- `removeUTRs_fromintrons` also returns the whole alignment as non-UTR (`aln[None:None]`),
- `sample_UTR3` raises `TypeError` on `None + 1` and is silently dropped as NotCompleted.

Such an alignment is counted twice, in two regions, and never as 3' UTR. Decide on one behaviour
(probably skip the alignment in all three apps).

## 3. Cosmetic and documentation drift

- `pipeline_filtering.md:6` says to set `DATA_HUMCHIMPORANG115`; `paths.py` and every script use
  `DATA_HUMCHIMPORANGOR114`. `pipeline_constraint_estimates.md` has the right name.
- `pipeline_filtering.md:21` says cds output is `filtered.fa`; scripts write `singlent_filtered.fa`.
- `bash_filtering3.sh` and `bash_submodel4.sh` echo `$REGION`, which is never set. `bash_get_humanseq.sh`
  and `bash_intron_filtering2.sh` echo "Running codon aligner".
- `libs.py:3` imports `define_app` from `scinexus.composable`; the `whole_genome_Orangutan/` copy uses
  `cogent3.app.composable`. It resolves in the Ensembl0.7.9 env but is an undocumented dependency not in
  `environment.yml`.
- `fittinglh_sm_4.py:14-20`: `singlentmodel_cds` and `singlentmodel_noncds` are identical. Singlent fits use
  `time_het="max"` with Orangutan as a discrete edge; trinuc fits use neither. Possibly deliberate.
- `mutation_rate_perregion.ipynb` is stale: references `paths.DATA_APES114` (no longer in `paths.py`), ran
  in a `delme` conda env, writes to a capitalised `Output_data/`, and duplicates the mutability model now in
  `estimating_mutrates.ipynb`. `pipeline_filtering.md:63` still points readers to it.
- `.gitignore` excludes every `*.csv`, so none of the output files or the duplicate problems are visible in
  git history.

## Suggested order of work
0. Section 0 first: re-download the five mis-masked regions; nothing downstream is meaningful until then.
1. Fix 1.1 (list reset) and 2.2 (alnstat swap), re-run `estimating_constraint_5.py`.
2. Re-run the `estimating_mutrates.ipynb` cell to regenerate a clean `mutability_motiflength3.csv` (1.2).
3. Settle the mean_sem CSV format (1.3) and the site-count definition (1.4), then re-run
   `plot_constraint.ipynb`, `plot_mutability.ipynb`, `Ud_estimation.ipynb`.
4. Fix 2.1 (CpG scale) before anything uses `get_scaled_lengths`.
5. Fix 2.3 and 2.4, then decide whether the human-sequence files need regenerating.
