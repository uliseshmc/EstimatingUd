#!/usr/bin/env python3
import cogent3
from cogent3 import get_app
import paths
import libs
import argparse
import os

SUBMODELS = ["trinuc", "singlent"]

def filtering_cds(substitutionmodel):
    loader = get_app("load_unaligned", moltype="dna")
    rename = libs.renamer_cds_unaligned()
    trim_stops = get_app("trim_stop_codons")
    codon_align = get_app("progressive_align", "codon", guide_tree="(Human:0.06,Chimpanzee:0.06,Orangutan:0.14)")
    if substitutionmodel == "singlent":
        omit_degs = get_app("omit_degenerates", moltype="dna", motif_length=1)
    elif substitutionmodel == "trinuc":
        omit_degs = get_app("omit_degenerates", moltype="dna", motif_length=3)
    else:
        raise ValueError("Trying to use a substitution model other than singlent, or trinuc")

    cds_app = loader + rename + trim_stops + codon_align + omit_degs
    return cds_app


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-submodel",
        "--substitutionmodel",
        type=str,
        required=True,
        choices=SUBMODELS,
        help="Substitution model (one of: " + ", ".join(SUBMODELS) + ")",
    )
    args = parser.parse_args()

    concat = get_app("concat", moltype="dna")    
    
    cds_app = filtering_cds(args.substitutionmodel)

    region = "cds/alldata_chrm22"
    folder_in = paths.DATA_HUMCHIMPORANGOR114 + region
    in_dstore = cogent3.open_data_store(folder_in, suffix='fa', mode='r')
    
    nonconcat_cds = [r for r in cds_app.as_completed(in_dstore[:], parallel=False) if r]
    cds_alns = concat(nonconcat_cds)

    region_out = "cds/chrm22"
    folder_out = paths.DATA_HUMCHIMPORANGOR114 + region_out
    os.makedirs(folder_out, exist_ok=True)

    if args.substitutionmodel == "singlent":
        label = "singlent_filtered"
    else:
        label = "trinucleotide_filtered"
        
    file_out = folder_out + "/" + label + ".fa"
    cds_alns.write(file_out)

    with open(folder_out + "/" + label + "_alnstat.txt", mode = "w") as out: 
        out.write("alignment length: " + str(len(cds_alns)))

if __name__ == "__main__":
    main()

