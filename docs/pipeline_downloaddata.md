# Installation and download

## The environment

Everything, downloads included, runs in one micromamba environment named `Udenv`
(Python 3.14, cogent3 2026.10.22, ensembl_tui 0.7.9). From the repo root:

```
micromamba create -f environment.yml
micromamba activate Udenv
pip install --no-deps -r requirements-udenv.txt
pip install --no-deps -e .
python -m ipykernel install --user --name Udenv
```

`requirements-udenv.txt` is the exact `pip freeze` of the working environment. It is
installed with `--no-deps` because two of its pins are ones pip's resolver would refuse:

- `scinexus==2026.9.12` although cogent3 2026.10.22 asks for 2026.10.2 or newer. The newer
  scinexus breaks `eti alignments` (`ensembl_tui/cli.py:633`, the app output becomes a
  `source_proxy`). Every scinexus name cogent3 imports exists in 2026.9.12 and the whole
  pipeline was checked under this pairing.
- `cogent3-h5seqs==0.9.1` although ensembl_tui 0.7.9 pins 0.7.3. Reading the existing
  `install/` (built with 0.7.3) works with 0.9.1. Writing a fresh install with 0.9.1 has
  not been verified; if `eti install` fails, fall back to
  `pip install "cogent3-h5seqs==0.7.3" "cogent3==2026.9.10"`.
- `duckdb<2`: duckdb 2.0 removed `Relation.fetchone()`, which `eti homologs` calls.

`pip check` reports both version conflicts. That is expected.

Tell the pipeline where the alignments are, once. Copy `.env.example` to `.env` in the repo
root and set the path there:

```
DATA_HUMCHIMPORANGOR114=/xdisk/masel/uliseshmc/EstimatingUd/HumChimpOranGor_10_114
```

`.env` is gitignored. Every `ud` command, notebook, `scripts/submit.sh` and SLURM wrapper reads
it when the `DATA_HUMCHIMPORANGOR114` variable is not already exported; exporting it in
`~/.bashrc` works too and takes precedence.

## Downloading and installing the alignments

```
micromamba activate Udenv
eti download -c HumChimpOran_10_114.cfg
cd HumChimpOran_10_114
eti install -d download -np <number of processors>
```

If a previous installation fails you can use `-f` to force the installation. `eti installed
-i install` shows the package versions an install was built with.

# Sampling the genomic regions

Each region below has a SLURM array wrapper in `scripts/download/` (one task per
chromosome, 1-22, X, Y). Submit them from the repo root:

```
scripts/submit.sh scripts/download/cds.sh
```

The wrappers activate `Udenv`, read `$DATA_HUMCHIMPORANGOR114`, and after each
`eti alignments` run print a mask report: how many fragments contain no `?` and what share
of the sampled human bases is masked. A download whose masking failed shows close to 100%
of fragments without `?` and close to 0% masked bases (this happened once for five regions;
see `docs/todo.md`). Check the `.out` logs before filtering.

## CDS

`eti homologs` exports one-to-one orthologous cds collections (unaligned; the pipeline
codon-aligns them in stage 1):

```
eti homologs -i install/ --outdir cds/alldata_chrm$chr --ref human --coord_names $chr
```

## Introns (masking cds and all repeats)

`eti alignments` samples the whole-genome alignment over protein-coding genes of the
reference species. It needs a file listing the biotypes to mask, one per line. List the
biotypes available for human with

```
eti species-summary -i install/ --species human
```

`cds_allAR_1column.txt` holds every unique entry of the first column of the repeats table
plus `cds` and `dust`.

<!-- Bug caution
eti species-summary outputs two columns for the repeats biotypes. We use all the unique entries in the first column.
-->

```
eti alignments -i install -od introns/alldata_chrm$chr --align_name "10_primates.epo" --ref human --mask cds_allAR_1column.txt --coord_names $chr
```

## Intronic ancestral repeats

Mask everything except ancestral repeats. We focus on LTRs, Type I Transposons/LINE,
Type I Transposons/SINE and Type II Transposons, listed in `focalAR_1column.txt`:

```
eti alignments -i install -od intronsAR/alldata_chrm$chr --align_name "10_primates.epo" --ref human --mask_shadow focalAR_1column.txt --coord_names $chr
```

## Intergenic regions

Intergenic sampling needs explicit coordinates (`--ref_coords`). Export the gene
metadata first:

```
eti dump-genes -i install --species human -od .
```

Then run `notebooks/01_location_inter_intragenic.ipynb` in `Udenv`. It reads
`homo_sapiens-114-gene_metadata.tsv`, takes the shadow of the gene spans and writes, per
chromosome, into `$DATA_HUMCHIMPORANGOR114/intergenic_coordinates/`:

- `chrom<chr>_intergenic.tsv` (all intergenic sequence, for intergenicAR)
- `chrom<chr>_proximal5_IG_coordinates.tsv`, `chrom<chr>_proximal3_IG_coordinates.tsv`,
  `chrom<chr>_distal_IG_coordinates.tsv`
- `chrom<chr>-intragenic.tsv` (gene spans, for reference)

Intergenic ancestral repeats (the neutral reference):

```
eti alignments -i install -od intergenicAR/alldata_chrm${chr} --align_name "10_primates.epo" --ref human --mask_shadow focalAR_1column.txt --ref_coords intergenic_coordinates/chrom${chr}_intergenic.tsv
```

Proximal 5', proximal 3' and distal intergenic regions mask all repeats with
`allAR_1column.txt` (same entries as `cds_allAR_1column.txt` without `cds` and `dust`):

```
eti alignments -i install -od proximal5IG/alldata_chrm${chr} --align_name "10_primates.epo" --ref human --mask allAR_1column.txt --ref_coords intergenic_coordinates/chrom${chr}_proximal5_IG_coordinates.tsv
eti alignments -i install -od proximal3IG/alldata_chrm${chr} --align_name "10_primates.epo" --ref human --mask allAR_1column.txt --ref_coords intergenic_coordinates/chrom${chr}_proximal3_IG_coordinates.tsv
eti alignments -i install -od distalIG/alldata_chrm${chr}    --align_name "10_primates.epo" --ref human --mask allAR_1column.txt --ref_coords intergenic_coordinates/chrom${chr}_distal_IG_coordinates.tsv
```

`--align_name` is matched as a glob against the installed alignment names; quote it or use
the exact name as above.

## Notes kept from earlier testing

### mafft alignments

To test alignment quality I used mafft to align the sequences. I started by removing gaps
and outputting the sequence collections to introns/mafft_aligned/unaligned_seqs, then

```
mafft --auto --treein tree.mafft sequence/.fa > ../sequence.fa
```

where the guide tree.mafft tree I used is

    1     2    0.06000    0.06000
    1     3    0.04000    0.10000

or in Newick format (3_Gorilla:0.1,(2_Chimpanzee:0.06,1_Human:0.06):0.04);

<!-- Bug caution

The eti alignments command was at one point unstable: it sometimes gave many warnings about
negative indexes, and in those cases the resulting files had no masked positions (no
question marks). The mask report printed by the download wrappers is there to catch this.

Code to debug this includes only focusing in one genomic region instead of multiples.
For this create a file "chrom22-selected.tsv" and manually indicate the coordinates needed.
My testing coordinates are 22:15915800-16141765. Then run

eti alignments -i install -od selected-subset21 --ref human --ref_coords chrom22-selected.tsv --mask_shadow focalAR_1column.txt --align_name "10_primates.epo"

-->

<!-- To check later

#Creates intragenic alignments for chromosome 22
eti alignments -i install -od test_intragenic --align_name "*primates*" --ref human --mask_shadow ancestralrepeats_list.txt --coord_names 22 --mask_ref --limit 10

#testing using the names of the alignment types instead of a file
eti alignments -i install --outdir test_intragenic --align_name "*primates*" --ref human --mask_shadow "SINE?/tRNA,SINE?,SINE/tRNA-Deu,SINE/tRNA,SINE/5S-Deu-L2,SINE/tRNA-RTE,SINE/MIR,SINE/Alu" --coord_names 22 --mask_ref --limit 2

#testing using the names of the alignment types instead of a file
eti alignments -i install --outdir test_intragenic --align_name "*primates*" --ref human --mask "cds" --coord_names 22 --mask_ref --limit 2

#Creates a sample of alignments of 20 genes from chromosome 1
eti alignments -i install/ --outdir aligns_demo/ --align_name '*primate*' --ref=human --limit=20 --coord_names=1

#Creates a sample of LINE alignments
eti alignments -i install --outdir aligns_line/chromosome1 --align_name '*primate*' --coord_names 1 --ref human --limit 10 --mask_shadow LINE --mask_ref

#Creates a sample of intron alignments
eti alignments -i install --outdir aligns_line/introns --align_name '*primate*' --coord_names 1 --ref human --limit 10 --mask cds --mask_ref

#Get gff3 annotations for humans
cd pathtodownloadannotations/
rsync -av rsync://ftp.ebi.ac.uk/ensemblorg/pub/release-112/gff3/homo_sapiens ./

-->
