#!/usr/bin/env python3
import cogent3
from cogent3 import get_app
import paths
import libs
import argparse
import os

CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y"]

def filtering_UTR5():
    loader = get_app("load_aligned", moltype="dna")
    rename_noncds = libs.renamer_noncds_aligned()
    get_UTR5 = libs.sample_UTR5()
    omit_gap_pos_app = get_app("omit_gap_pos", moltype="dna")
    omit_degs_noncds = get_app("omit_degenerates", moltype="dna", motif_length=3)
    UTR5_app = loader + rename_noncds + get_UTR5 + omit_gap_pos_app + omit_degs_noncds

    return UTR5_app

def filtering_UTR3():
    loader = get_app("load_aligned", moltype="dna")
    rename_noncds = libs.renamer_noncds_aligned()
    get_UTR3 = libs.sample_UTR3()
    omit_gap_pos_app = get_app("omit_gap_pos", moltype="dna")
    omit_degs_noncds = get_app("omit_degenerates", moltype="dna", motif_length=3)
    UTR3_app = loader + rename_noncds + get_UTR3 + omit_gap_pos_app + omit_degs_noncds

    return UTR3_app 

def filtering_nonUTR():
    loader = get_app("load_aligned", moltype="dna")
    rename_noncds = libs.renamer_noncds_aligned()
    remove_UTR = libs.removeUTRs_fromintrons()
    omit_gap_pos_app = get_app("omit_gap_pos", moltype="dna")
    omit_degs_noncds = get_app("omit_degenerates", moltype="dna", motif_length=3)
    nonUTR_app = loader + rename_noncds + remove_UTR + omit_gap_pos_app + omit_degs_noncds

    return nonUTR_app

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-chrm",
        "--chromosome",
        type=str,
        required=True,
        choices=CHROMOSOMES,
        help="Chromosome to process (one of: " + ", ".join(CHROMOSOMES) + ")",
    )
    args = parser.parse_args()
    concat = get_app("concat", moltype="dna")

    ## Process UTR5
    UTR5_app = filtering_UTR5()
    
    region = "introns/alldata_chrm" + args.chromosome
    folder_in = paths.DATA_HUMCHIMPORANG115 + region
    in_dstore = cogent3.open_data_store(folder_in, suffix='fa', mode='r')
    
    nonconcat_introns = [r for r in UTR5_app.as_completed(in_dstore[:], parallel=False) if r]
    introns_alns = concat(nonconcat_introns)

    region_out = "introns5UTR/chrm" + args.chromosome
    folder_out = paths.DATA_HUMCHIMPORANG115 + region_out
    os.makedirs(folder_out, exist_ok=True)

    file_out = folder_out + "/trinucleotide_filtered.fa"
    introns_alns.write(file_out)
    #print("Processed region:", region)

    ## Process UTR3
    UTR3_app = filtering_UTR3()
    
    nonconcat_introns = [r for r in UTR3_app.as_completed(in_dstore[:], parallel=False) if r]
    introns_alns = concat(nonconcat_introns)

    region_out = "introns3UTR/chrm" + args.chromosome
    folder_out = paths.DATA_HUMCHIMPORANG115 + region_out
    os.makedirs(folder_out, exist_ok=True)

    file_out = folder_out + "/trinucleotide_filtered.fa"
    introns_alns.write(file_out)
    #print("Processed region:", region)

    ## Process non-UTR introns
    nonUTR_app = filtering_nonUTR()
    
    nonconcat_introns = [r for r in nonUTR_app.as_completed(in_dstore[:], parallel=False) if r]
    introns_alns = concat(nonconcat_introns)
    
    region_out = "introns_nonUTR/chrm" + args.chromosome
    folder_out = paths.DATA_HUMCHIMPORANG115 + region_out
    os.makedirs(folder_out, exist_ok=True)

    file_out = folder_out + "/trinucleotide_filtered.fa"
    introns_alns.write(file_out)
    #print("Processed region:", region)

if __name__ == "__main__":
    main()

