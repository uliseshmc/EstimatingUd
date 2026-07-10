#!/usr/bin/env python3
import cogent3
from cogent3 import get_app
import paths
import libs
import argparse

CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y"]

def filtering_cds():
    loader = get_app("load_unaligned", moltype="dna")
    rename = libs.renamer_cds_unaligned()
    trim_stops = get_app("trim_stop_codons")
    codon_align = get_app("progressive_align", "codon", guide_tree="(Human:0.06,Chimpanzee:0.06,Orangutan:0.14)")
    cds_app = loader + rename + trim_stops + codon_align

    return cds_app

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

    cds_app = filtering_cds()

    region = "cds/chrm" + args.chromosome
    folder_in = paths.DATA_HUMCHIMPORANG115 + region
    in_dstore = cogent3.open_data_store(folder_in, suffix='fa', mode='r')

    nonconcat_cds = [r for r in cds_app.as_completed(in_dstore[:], parallel=False) if r]
    cds_alns = concat(nonconcat_cds)

    file_out = paths.DATA_HUMCHIMPORANG115 + region + "/filtered.fa"
    cds_alns.write(file_out)
    #print("Processed region:", region)

if __name__ == "__main__":
    main()

