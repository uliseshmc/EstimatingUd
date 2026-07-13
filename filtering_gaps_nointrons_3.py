#!/usr/bin/env python3
import cogent3
from cogent3 import get_app
import paths
import libs
import argparse
import os


REGIONS = ["intergenicAR", "intronsAR", "distalIG", "proximal5IG", "proximal3IG"]
CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y"]

def filtering_noncds():
    loader = get_app("load_aligned", moltype="dna")
    rename_noncds = libs.renamer_noncds_aligned()
    omit_gap_pos_app = get_app("omit_gap_pos", moltype="dna")
    noncds_app = loader + rename_noncds + omit_gap_pos_app

    return noncds_app

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-reg",
        "--region",
        type=str,
        required=True,
        choices=REGIONS,
        help="Region to process (one of: " + ", ".join(REGIONS) + ")",
    )
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
    
    region = args.region + "/alldata_chrm" + args.chromosome
    folder_in = paths.DATA_HUMCHIMPORANG115 + region
    in_dstore = cogent3.open_data_store(folder_in, suffix='fa', mode='r')
    
    nonconcat_noncds = [r for r in noncds_app.as_completed(in_dstore[:], parallel=False) if r]
    noncds_alns = concat(nonconcat_noncds)

    region_out = args.region + "/chrm" + args.chromosome
    folder_out = paths.DATA_HUMCHIMPORANG115 + region_out
    os.makedirs(folder_out, exist_ok=True)
    
    file_out = folder_out + "/filtered.fa"
    noncds_alns.write(file_out)
    #print("Processed region:", region)

if __name__ == "__main__":
    main()

