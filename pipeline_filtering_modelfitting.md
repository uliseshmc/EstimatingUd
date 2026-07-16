# Preparation

This instructions assume that pipeline_downloaddata.md is finished. 
Instructions for installation of the conda package UdChimpHumOran can be found there.

Make sure that the variable DATA_HUMCHIMPORANG115 in the file paths.py points at the location where you downloaded the alignments. 

## Filtering cds

eti homologs (described in pipeline_download.md) creates a sequence collection of cds orthologous cds regions across the Chimps, Orangutans and Human. To align these sequences at the codon level we use

```
conda activate UdChimpHumOran

python3 codon_aligner_1.py -chrm $chr
```

This script will align and remove stop codons from the cds sequences. It will output the conncatenetad whole chromosome alignment in the file filtered.fa. This and all other filtered files are ocated in the same folder as its corresponding chromosome datastore

For trinucleotide models, trinucletide sites with any gap have to be filtered out from the alignment. This is done by the script 

```
conda activate UdChimpHumOran

python3 codon_aligner_trinucs_1.py -chrm $chr
```
where $chr is the chromosome stableid (1,2,..22, X, Y)

This and all trinucleotide filtering are output in a file named trinucleotide_filtered.fa.

## Filtering gapped sites (in regions other than cds and introns)

The alignment of Chimps, Orangutans and Human contain many gap sites corresponding to data in any of the other 7 aligned species. Also, the data is distributed accross many files that take up more than 20 TB of storage. To solve this issue, I filtered out all gaped sites and store the whole chromosome alignment into a single concatenated alignment in the file filtered.fa. This is done by.

```
conda activate UdChimpHumOran

python3 filtering_gaps_nointrons_3.py -reg <insert $region -chrm $chr
```

where available regions are "intergenicAR", "intronsAR", "distalIG", "proximal5IG", or "proximal3IG".

Trinucleotide filtering is achieved by

```
conda activate UdChimpHumOran

python3 filtering_trinucs_gaps_nointrons_3.py -reg $region -chrm $chr
```


## Filtering introns

I divided intron sequences into 5'UTR, 3'UTR and nonUTR regions. Then I filtered out gapped sites (See section above). This is done by 

```
conda activate UdChimpHumOran

python3 filtering_gaps_introns_2.py -chrm $chr
```

For trinucleotide models, use

```
conda activate UdChimpHumOran

python3 filtering_trinucs_gaps_introns_2.py -chrm $chr
```

## Bash mode

If using a cluster that runs under a SLURM system, you can use bash*.sh files to run filtering for all regions through all seqids.

Sometimes SLURM fails. Please double checked the .err files to check for common errors. A typical flagged error is AttributeError: 'NotCompleted' object has no attribute 'write'. This usually happens because the length of the alignment is 0. 
To make a list of all other errors I used

```
find . -maxdepth 1 -type f -name "*.err" -size +0c -exec grep -FL "AttributeError: 'NotCompleted' object has no attribute 'write'" {} \; > runs_with_errors.txt
```
This will output a list of all .err files into runs_with_errors.txt.

# Dealing with storage limitations

The eti command originally creates a file for each contig. This result in huge storage requirements when downloading the whole genome. The data is divided into two folder categories.

$region/alldata_chrm$chr

stores the outputs of the eti command. 

$region/chrm$chr

stores the filtered data resulting from codon_aligner_1.py, filtering_gaps_introns_2.py and filtering_gaps_nointrons.py (and their respective trinucleotide versions). 

If running into data storage limitations, I recommend to zip all data using 

```
tar -cf - "$region" | pigz -p $#processors > "$region.tar.gz"
```

and then delete all the $region/alldata_chrm$chr folders.