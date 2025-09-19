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

```
pip install ensembl_tui
pip install "cogent3==2025.7.10a5" "cogent3_h5seqs==0.5.0" -U
```
Please check that the package versions are cogent3==2025.7.10a5, cogent3-h5seqs==0.5.0 and ensembl_tui==0.4.3 using
```
pip freeze
```


## Installation of other libraries
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
-np is number of process

if a previous installation fails you can use -f to force the installation

# Analysis of genomic regions

Currently sampling only chromosome 22 to speed up analysis

## CDS

We start by sampling homologs alignments

```
eti homologs -i install/ --outdir aligns_cds --ref human --coord_names 22 
```

The resulting alignments were made using synteny so we know that specie's regions are orthologous.

However, the alignment has not been done at the codon level and thus it suffers from alignment errors.

We first treat the data with a codon aligner using the python notebook codon_aligner.ipynb



## Intergenic regions (not ancestral repeats)

First we need to create a file with genomic coordinates 
```
eti dump-genes -i install --species human -od .
```

Then run the python notebook "location_intergenic_regions.ipynb" 

This program uses the file produced by dump-genes to extract the location coordinates of intergenic regions

It produces the file "chrom22-intergenic.tsv". We use this file to indicate where eti has to look for alignments

The second part that eti requires is a file with the biotypes we want to do a mask_shadow. The next command shows the biotypes available

```
eti species-summary -i install/ --species human
```

<!-- Bug caution
eti spcies-summary outputs two columns for the reapeats bitype. We use all the unique entries in the first column.
We could also work with the second column but it masks less sequences than the first, so I'm using the first column.
-->

We want to mask in our alignments all the data coming from reapeats sequences. 

For this, we make a file allAR_1column.txt containing all the unique biotypes in the first column of table repeats from the eti species-summary command. Each biotype per line.

Then we sample alignments with ancestral repeat sequences masked using

```
eti alignments -i install -od intergenic --align_name 10_primates* --ref human --mask allAR_1column.txt --ref_coords chrom22-intergenic.tsv
```

## Intergenic ancestral repeats


We will focus in ancestral repeats that belongs to LTRs, Type I Transposons/LINE, Type I Transposons/SINE and Type II Transposons. The file ancestralrepeats_1column.txt containes these bitypes.



Then we extract alignments with everything masked but these biotypes using
```
eti alignments -i install -od test_intergenicAR_1column --align_name 10_primates* --ref human --mask_shadow ancestralrepeats_1column.txt --ref_coords chrom22-intergenic.tsv
```



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
