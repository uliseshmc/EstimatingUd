# Filtering

These instructions assume `docs/pipeline_downloaddata.md` is finished: `Udenv` is
installed with the package (`pip install -e .`), `DATA_HUMCHIMPORANGOR114` points at the
alignment root, and the raw `eti` output sits under `<region>/alldata_chrm<chr>/`.

```
micromamba activate Udenv
export DATA_HUMCHIMPORANGOR114=/path/to/HumChimpOranGor_10_114
```

Every command takes `--chrm` (one of `1`-`22`, `X`, `Y`) and `--submodel`
(`singlent` or `trinuc`). The submodel only changes the motif length passed to
`omit_degenerates` (1 or 3) and the output filename: `singlent_filtered.fa` or
`trinucleotide_filtered.fa`, each with a `_alnstat.txt` holding the alignment length.
Outputs land in `<region>/chrm<chr>/` next to the region's raw data store.

## Stage 1: cds

`eti homologs` produces unaligned collections of orthologous cds. Stage 1 renames the
species, trims stop codons, codon-aligns with the guide tree
`(Human:0.06,Chimpanzee:0.06,Orangutan:0.14)`, drops degenerate columns (masked `?`,
gaps, non-nucleotide symbols) and concatenates the genes of a chromosome:

```
ud align-cds --chrm $chr --submodel singlent
ud align-cds --chrm $chr --submodel trinuc
```

For the trinucleotide models, every trinucleotide with a degenerate position is removed
rather than every column.

## Stage 2: introns

The intron alignments are split into introns5UTR (from the transcript start to the first
masked position in the Human sequence), introns3UTR (after the last masked position) and
introns_nonUTR (between the two). All three read `introns/alldata_chrm<chr>` and write to
their own region folder:

```
ud filter-introns --chrm $chr --submodel singlent
ud filter-introns --chrm $chr --submodel trinuc
```

## Stage 3: the other non-cds regions

intergenicAR, intronsAR, distalIG, proximal5IG and proximal3IG only need renaming and
degenerate filtering:

```
ud filter-noncds --chrm $chr --submodel singlent
ud filter-noncds --chrm $chr --submodel trinuc
```

## Human sequence per region

The mutability model (`notebooks/03_estimating_mutrates.ipynb`) needs the unaligned
Human sequence of every region, concatenated per chromosome and truncated to whole
trinucleotides. The intron regions are sampled exactly as in stage 2 first.

```
ud human-seq --chrm $chr
```

This writes `<region>/chrm<chr>/human_seq.fa` for all nine regions.

## SLURM

`scripts/slurm/` holds one `--array=1-24` wrapper per stage (task ID to chromosome),
running the singlent then the trinuc invocation. Submit from the repo root so that
`SLURM_SUBMIT_DIR` is the repo and `logs/` is created:

```
scripts/submit.sh scripts/slurm/01_align_cds.sh
scripts/submit.sh scripts/slurm/02_filter_introns.sh
scripts/submit.sh scripts/slurm/03_filter_noncds.sh
scripts/submit.sh scripts/slurm/03b_human_seq.sh
```

Logs are `logs/<job name>_<array job>_<task>.out|err`. When a task fails, read the `.err`.
`AttributeError: 'NotCompleted' object has no attribute 'write'` means every fragment of a
region was rejected (alignment of length 0) and can usually be ignored. List the other
failures with

```
find logs -maxdepth 1 -type f -name "*.err" -size +0c -exec grep -FL "AttributeError: 'NotCompleted' object has no attribute 'write'" {} \; > runs_with_errors.txt
```

A `FileNotFoundError` for `cds/chrmY` is expected: there is no cds alignment for Y.

## Dealing with storage limitations

Once stages 1 to 3 and `ud human-seq` have run for a region, its raw
`alldata_chrm<chr>` folders can be archived and deleted:

```
tar -cf - "$region" | pigz -p <processors> > "$region.tar.gz"
```
