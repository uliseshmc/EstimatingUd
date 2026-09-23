# Preparation

This instructions assume that pipeline_downloaddata.md is finished. 
Instructions for installation of the conda package Ensembl0.7.9 can be found there.

Make sure that the variable DATA_HUMCHIMPORANG115 in the file paths.py points at the location where you downloaded the alignments. 

Start by activating the environment
```
conda activate Ensembl0.7.9
```

## Filtering cds

eti homologs (described in pipeline_download.md) creates a sequence collection of cds orthologous cds regions across the Chimps, Orangutans and Human. To align these sequences at the codon level we use

```
python3 codon_aligner_1.py -chrm $chr -submodel singlent
```

This script will align and remove stop codons from the cds sequences. It will also remove columns with degenerate characters (maksed: "?", gaps: "-", non neucleotide symbols). It will output the conncatenetad whole chromosome alignment in the file filtered.fa. This and all other filtered files are located in the same folder as its corresponding chromosome datastore

For trinucleotide models, trinucletide sites with any gap have to be filtered out from the alignment. This is done by the script 

```
python3 codon_aligner_1.py -chrm $chr -submodel trinuc
```
where $chr is the chromosome stableid (1,2,..22, X, Y)

This and all trinucleotide filtering are output in a file named trinucleotide_filtered.fa.

## Filtering gapped sites (in regions other than cds and introns)

For the sinle nucleotide model

```
python3 filtering_gaps_nointrons_3.py -chrm $chr -submodel singlent
```

For the trinucleotide models

```
python3 filtering_gaps_nointrons_3.py -chrm $chr -submodel trinuc
```


## Filtering introns

I divided intron sequences into 5'UTR, 3'UTR and nonUTR regions. Then I filtered out gapped sites (See section above). This is done by 

```
python3 filtering_gaps_introns_2.py -chrm $chr -submodel singlent
```

For trinucleotide models, use

```
python3 filtering_gaps_introns_2.py -chrm $chr -submodel trinuc
```

## Taking only human sequence

Our mutation model (see mutation_rate_perregion.ipynb) considers genomic regions to account for mutation heterogeneity.
To do this. we need to get the unconncatenated human sequence. This is done by

```
python3 get_human_sequence.py -chrm $chr -mutmotif 3
```

We use the Oman et al model which accounts for trinuleotide context. Thus we set mutmotif to 3

## Bash mode

If using a cluster that runs under a SLURM system, you can use bash*.sh files to run filtering for all regions through all seqids.

Sometimes SLURM fails. Please double checked the .err files to check for common errors. A typical flagged error is AttributeError: 'NotCompleted' object has no attribute 'write'. This usually happens because the length of the alignment is 0. 
To make a list of all other errors I used

```
find . -maxdepth 1 -type f -name "*.err" -size +0c -exec grep -FL "AttributeError: 'NotCompleted' object has no attribute 'write'" {} \; > runs_with_errors.txt
```
This will output a list of all .err files into runs_with_errors.txt.

# Dealing with storage limitations

If running into data storage limitations, I recommend running the filtering and then to zip alldata folders using 

```
tar -cf - "$region" | pigz -p $#processors > "$region.tar.gz"
```

and then delete all the $region/alldata_chrm$chr folders.