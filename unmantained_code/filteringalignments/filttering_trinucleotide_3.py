#!/usr/bin/env python3
import cogent3
from cogent3 import get_app
from EstimatingUd import paths
from EstimatingUd import libs
import argparse

def filtering_cds():
    loader_cds = get_app("load_aligned", moltype="dna")
    omit_gap_pos_app = get_app("omit_gap_pos", moltype="dna")
    omit_degs_cds = get_app("omit_degenerates", moltype="dna", motif_length=3)
    cds_process = loader_cds + omit_gap_pos_app + omit_degs_cds

    return cds_process


def filtering_noncds():
    loader = get_app("load_aligned", moltype="dna")
    rename_noncds = libs.renamer_noncds_aligned()
    omit_gap_pos_app = get_app("omit_gap_pos", moltype="dna")
    omit_degs_noncds = get_app("omit_degenerates", moltype="dna", motif_length=3)
    noncds_app = loader + rename_noncds + omit_gap_pos_app + omit_degs_noncds

    return noncds_app

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-chrm", "--chromosome", type=str, required=True, help='Chromosome to process')
    args = parser.parse_args()
    
    concat = get_app("concat", moltype="dna")

    folder_in = paths.DATA_HUMCHIMPORANG115 + 'cds/chrm' + args.chromosome + '/codon_aligned'
    in_dstore = cogent3.open_data_store(folder_in, suffix='fa', mode='r')

    cds_app = filtering_cds()
    nonconcat_cds = [r for r in cds_app.as_completed(in_dstore[:], parallel=False) if r]
    cds_alns = concat(nonconcat_cds)

    file_out = paths.DATA_HUMCHIMPORANG115 + 'cds/chrm' + args.chromosome + '/codon_aligned/trinucleotide_filtered.fa'
    cds_alns.write(file_out)

    #run for  chromosome 22 finished 3UTR
    noncdsregions = ["intergenicAR/chrm" + args.chromosome, 
                     "introns/chrm" + args.chromosome + "/noUTRs", 
                     "introns/chrm" + args.chromosome + "/5UTR", 
                     "introns/chrm" + args.chromosome + "/3UTR", 
                     "intronsAR/chrm" + args.chromosome, 
                     "distalIG/chrm" + args.chromosome,
                     "proximal5IG/chrm" + args.chromosome, 
                     "proximal3IG/chrm" + args.chromosome]

    noncds_app = filtering_noncds()
    for region in noncdsregions:
        folder_in = paths.DATA_HUMCHIMPORANG115 + region
        in_dstore = cogent3.open_data_store(folder_in, suffix='fa', mode='r')
        
        nonconcat_noncds = [r for r in noncds_app.as_completed(in_dstore[:], parallel=False) if r]
        noncds_alns = concat(nonconcat_noncds)

        file_out = paths.DATA_HUMCHIMPORANG115 + region + "/trinucleotide_filtered.fa"
        noncds_alns.write(file_out)


if __name__ == "__main__":
    main()
