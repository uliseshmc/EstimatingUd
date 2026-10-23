#!/usr/bin/env python3
import cogent3
from cogent3 import get_app
import paths
import libs
import argparse
import os

CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y"]
REGIONS = ["intergenicAR", "intronsAR", "distalIG"]

def filtering_noncds():
    loader = get_app("load_aligned", moltype="dna")
    rename_noncds = libs.renamer_noncds_aligned()
    omit_degs = get_app("omit_degenerates", gap_is_degen=False, moltype="dna", motif_length=1)
    remove_gaps = get_app("omit_gap_pos", allowed_frac=0.90)
    
    noncds_app = loader + rename_noncds + omit_degs + remove_gaps

    return noncds_app

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

    noncds_app = filtering_noncds()
    concat = get_app("concat", moltype="dna")

    for genomic_region in REGIONS:
        relative_folder_in = genomic_region + "/alldata_chrm" + args.chromosome
        folder_in = paths.DATA_HUMCHIMPORANGOR114 + relative_folder_in
        in_dstore = cogent3.open_data_store(folder_in, suffix='fa', mode='r')
        
        nonconcat_noncds = [r for r in noncds_app.as_completed(in_dstore[:], parallel=False) if r]
        noncds_alns = concat(nonconcat_noncds)

        relative_folder_out = genomic_region + "/chrm" + args.chromosome
        folder_out = paths.DATA_HUMCHIMPORANGOR114 + relative_folder_out
        os.makedirs(folder_out, exist_ok=True)
        
        label = "singlent_nomasked"
        
        file_out = folder_out + "/" + label + ".fa"
        noncds_alns.write(file_out)
        
        with open(folder_out + "/" + label + "_alnstat.txt", mode = "w") as out: 
            out.write("alignment length: " + str(len(noncds_alns)))

if __name__ == "__main__":
    main()

