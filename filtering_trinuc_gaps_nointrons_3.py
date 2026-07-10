#!/usr/bin/env python3
import cogent3
from cogent3 import get_app
import paths
import libs
import argparse

REGIONS = ["intergenicAR", "intronsAR", "distalIG", "proximal5IG", "proximal3IG"]
CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y"]

#Same as filtering_gaps_nointrons_3.py but omit_degs_noncds is added to considered only aligned trinucleotides with no gaps.
def filtering_noncds():
    loader = get_app("load_aligned", moltype="dna")
    rename_noncds = libs.renamer_noncds_aligned()
    omit_gap_pos_app = get_app("omit_gap_pos", moltype="dna")
    omit_degs_noncds = get_app("omit_degenerates", moltype="dna", motif_length=3)
    noncds_app = loader + rename_noncds + omit_gap_pos_app + omit_degs_noncds

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
    
    region = args.region + "/chrm" + args.chromosome
    folder_in = paths.DATA_HUMCHIMPORANG115 + region
    in_dstore = cogent3.open_data_store(folder_in, suffix='fa', mode='r')
    
    nonconcat_noncds = [r for r in noncds_app.as_completed(in_dstore[:], parallel=False) if r]
    noncds_alns = concat(nonconcat_noncds)

    file_out = paths.DATA_HUMCHIMPORANG115 + region + "/trinucleotide_filtered.fa"
    noncds_alns.write(file_out)
    #print("Processed region:", region)

if __name__ == "__main__":
    main()

