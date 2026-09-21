# Preparation

This instructions assume that pipeline_downloaddata.md is finished. 
There you can also find instructions for installation of the conda package Ensembl0.7.7.

Make sure that the variable DATA_HUMCHIMPORANGOR114 in the file paths.py points at the location where you downloaded the alignments. 

Start by activating the environment
```
conda activate Ensembl0.7.7
```

# Filtering eti alignment output

## CDS

eti homologs (described in pipeline_download.md) creates a sequence collection of cds orthologous cds regions across the Chimps, Orangutans and Human. To align these sequences at the codon level we use

```
python3 codon_aligner_1.py -submodel singlent
```

This script will align, remove stop codons, and remove degenerate characters from the cds sequences. 

For trinucleotide models, trinucletide sites with any degenerate character have to be filtered out from the alignment. This is done by the script 

```
python3 codon_aligner_1.py -submodel trinuc
```

These and all other scripts will output the conncatenetad whole chromosome alignment in the file singlent_filtered.fa or trinucleotide_filtered.fa. These scripts also output *_alnsst.txt file reporting the alignment lenght of each file. All files are located in the same folder as its corresponding chromosome datastore.

## Introns

I divided intron sequences into 5'UTR, 3'UTR and nonUTR regions. Then I filtered out degenerate sites (See section above). This is done by 

```
python3 filtering_gaps_introns_2.py -submodel singlent
```

For trinucleotide models, use

```
python3 filtering_gaps_introns_2.py -submodel trinuc
```

## Other regions

For all other regions

```
python3 filtering_gaps_nointrons_3.py -reg $region -submodel singlent
```

where available regions are "intergenicAR", "intronsAR", "distalIG", "proximal5IG", or "proximal3IG".

Trinucleotide filtering is achieved by

```
python3 filtering_gaps_nointrons_3.py -reg $region -submodel trinuc
```

# Bash mode

If using a cluster that runs under a SLURM system, you can use bash*.sh files to run filtering for all regions through all seqids.

Sometimes SLURM fails. Please double checked the .err files to check for common errors. A typical flagged error is AttributeError: 'NotCompleted' object has no attribute 'write'. This usually happens because the length of the alignment is 0. 
To make a list of all other errors I used

```
find . -maxdepth 1 -type f -name "*.err" -size +0c -exec grep -FL "AttributeError: 'NotCompleted' object has no attribute 'write'" {} \; > runs_with_errors.txt
```
This will output a list of all .err files into runs_with_errors.txt.

# Dealing with storage limitations

The eti alignment output files for the whole genome take up more than 20 TB of storage. This is because the alignment of Chimps, Orangutans and Human contained also 7 other primate species. This results in many gap sites. Also by masking sites we duplicate positions. For example intorns and intronsAR have the same number of sites, one being the shadow of the other. The eti command also creates a file for each contig making storage requirments bigger. 

To release storage you can run one eti alignment at a time, then filtered the output and then delete the tar the alldata_chrm$chrm folders. That is:

The data is divided into two folder categories.

$region/alldata_chrm$chr

stores the outputs of the eti command. 

$region/chrm$chr

stores the filtered data resulting from codon_aligner_1.py, filtering_gaps_introns_2.py and filtering_gaps_nointrons.py (and their respective trinucleotide versions). 

If running into data storage limitations, I recommend to zip all data using 

```
tar -cf - "$region" | pigz -p $#processors > "$region.tar.gz"
```

and then delete all the $region/alldata_chrm$chr folders.