# Installation
## Conda Environment setup
```
conda create -n delme python=3.13
conda activate delme
```

## Installing cogent3 and ensembl_tui

<!-- I'm not doing any of this now
through github
```
mkdir -p ~/repos
cd ~/repos
cd ../
git clone git@github.com:cogent3/ensembl_tui.git
cd ensembl_tui
pip install -e ".[dev]"
git clone git@github.com:cogent3/cogent3.git
cd cogent3
pip install -e ".[dev]"
```

through flit
```
pip install flit
install -s --python `which python
```
-->

We start by installing ensembl tui and cogent3 packagest
```
pip install ensembl_tui
pip install "cogent3==2025.7.10a5" "cogent3_h5seqs==0.5.0" -U
```
Please check that the package versions are cogent3==2025.7.10a5, cogent3-h5seqs==0.5.0 and ensembl_tui==0.4.3 using
```
pip freeze
```


## Installation of other libraries
We also need this libraries
```
pip install pandas
pip install phyilim
```

## Downloading and installing the alignments
```
conda activate delme
eti download -c primates10_114.cfg
cd primates10_114
eti install -d download -np <insert number>
```
-np is number of processors

if a previous installation fails you can use -f to force the installation

# Analysis of genomic regions

Currently sampling only chromosome 22 to speed up analysis. In the final analysis I should do all analysis for the whole genome data

## CDS

We start by sampling homologs alignments

```
eti homologs -i install/ --outdir cds --ref human --coord_names 22 
```

The resulting alignments were made using synteny so we know those regions are orthologous.

However, the alignment has not been done at the codon level and thus it suffers from alignment errors. We make a directory inside the cds folder to store codon aligned sequences

```
mkdir  cds/codon_aligned
```

Now we treat the data with a codon aligner using the python notebook codon_aligner.ipynb

Then we fit a substitution model to the new aligned sequences using CDS_substitutionmodel.ipynb

## Introns (masking ancestral repeats)

Here we use the eti alignments command. This command looks by default genomic coordinates corresponding to protein coding genes for the reference species. 

```
eti alignments -i install -od introns --align_name 10_primates* --ref human --mask cds_allAR_1column.txt --coord_names 22 
```
Here, we mask cds and and all ancestral repeatsto by using the file cds_allAR_1column.txt

Then we fit a substitution model to the new aligned sequences using introns_substitutionmodel.ipynb

## Intronic ancestral repeats

We again use the eti alignments command. Instead of masking ancestral repeats now we using mask_shadow to mask everything but the ancestral repeats.

```
eti alignments -i install -od intronsAR --align_name 10_primates* --ref human --mask_shadow ancestralrepeats_1column.txt --coord_names 22 
```

## Intragenic ancestral repeats

We use again the eti alignments command. However, this time we want our query to be inside of the gene. Thus we need to provide the genomic coordinates using the --ref_coords flag.

First we need to know the intragenomic coordinates. We use the command eti dump-genes to create a file with such coordinates homo_sapiens-114-gene_metadata.tsv

```
eti dump-genes -i install --species human -od .
```

The file homo_sapiens-114-gene_metadata.tsv contain too many columns not needed in upstream analysis. For this we run "location_inter_intragenic.ipynb" to simplify the data and output it on chrom22-intragenic.tsv

Second, we need a list of ancestral elements that we want to analyse. We will focus on ancestral repeats of type LTRs, Type I Transposons/LINE, Type I Transposons/SINE and Type II Transposons. We write these categories on the file ancestralrepeats_1column.txt one per line.

We use the command eti alignments to produce alignments of the intergenic regions. We will use the flag  use the the flag --mask_shadow to mask everything in our alignments but the ancestral repeats types included in ancestralrepeats_1column.txt. We will also use the flag --ref_coords  to indicate that we want to sample only from the intragenomic regions contained at chrom22-intergenic.tsv

```
eti alignments -i install -od intragenicAR --align_name 10_primates* --ref human --mask_shadow ancestralrepeats_1column.txt --ref_coords chrom22-intragenic.tsv
```


## Intergenic regions (masking ancestral repeats)

We need to know the intergenetic coordinates. For this we run the python notebook "location_inter_intragenic.ipynb". This notebook takes the file homo_sapiens-114-gene_metadata.tsv created by the eti-dump command and outputs a file with intergenic coordinates chrom22-intergenic.tsv. We use this file to indicate eti where to look for the alignments.

The second file that eti requires is one that indicates the biotypes we want to mask. To find all the biotypes available use

```
eti species-summary -i install/ --species human
```

<!-- Bug caution
eti spcies-summary outputs two columns for the reapeats bitypes. We use all the unique entries in the first column.
We could also work with the second column but it masks less sequences than the first, so I'm using the first column.
-->

We want to mask all the data coming from reapeats sequences. 

For this, we make a file 'allAR_1column.txt' containing all the unique biotypes in the first column of table repeats from the eti species-summary command. Each biotype per line.

Then we sample alignments with ancestral repeat sequences masked using

```
eti alignments -i install -od intergenic --align_name 10_primates* --ref human --mask allAR_1column.txt --ref_coords chrom22-intergenic.tsv
```

We then fit a substitution model using CDS_substitutionmodel.ipynb

## Intergenic ancestral repeats

Here we use the same steps as in the section "Intergenic ancestral repeats" but instead we use chrom22-intergenic.tsv to indicate that we want to sample from intergenic regions

```
eti alignments -i install -od intergenicAR --align_name 10_primates* --ref human --mask_shadow ancestralrepeats_1column.txt --ref_coords chrom22-intergenic.tsv
```
We then fit a substitution model using intergenicAR_substitutionmodel.ipynb


<!-- Bug caution

The previous command is quite unstable. Sometimes it gives a lot of warnings stating that the user is attempting to use negative indexes. In these cases the resulting files inside of test_intergenic_1column do not have any masked positions (No question marks)

Code to debug this includes only focusing in one genomic region instead of multiples.
For this create a file "chrom22-selected.tsv" and manually indicate the coordinates needed

My testing coordinates are 22:15915800-16141765


Then run
```
eti alignments -i install -od selected-subset21 --ref human --ref_coords chrom22-selected.tsv --mask_shadow ancestralrepeats_1column.txt --align_name 10_primates*
```

to generate the alignment

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
