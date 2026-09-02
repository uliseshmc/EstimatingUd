#!/usr/bin/env python3
import cogent3
from cogent3 import get_app
import whole_genome_Orangutan.paths as paths
import whole_genome_Orangutan.libs as libs
import argparse
import os

SUBMODELS = ["trinuc", "singlent"]
REGIONS = ["intergenicAR", "intronsAR", "distalIG", "proximal5IG", "proximal3IG"]

def filtering_noncds(trinucleotide):
    loader = get_app("load_aligned", moltype="dna")
    rename_noncds = libs.renamer_noncds_aligned()
    if trinucleotide == "singlent":
        omit_degs = get_app("omit_degenerates", moltype="dna", motif_length=1)
    elif trinucleotide == "trinuc":
        omit_degs = get_app("omit_degenerates", moltype="dna", motif_length=3)
    else:
        raise ValueError("Trying to use a substitution model other than singlent, or trinuc")
    noncds_app = loader + rename_noncds + omit_degs

    return noncds_app

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

    noncds_app = filtering_noncds(args.substitutionmodel)
    concat = get_app("concat", moltype="dna")

    for genomic_region in REGIONS:
        region = genomic_region + "/alldata_chrm22"
        folder_in = paths.DATA_HUMCHIMPORANG115 + region
        in_dstore = cogent3.open_data_store(folder_in, suffix='fa', mode='r')
        
        nonconcat_noncds = [r for r in noncds_app.as_completed(in_dstore[:], parallel=False) if r]
        noncds_alns = concat(nonconcat_noncds)

        region_out = genomic_region + "/chrm22"
        folder_out = paths.DATA_HUMCHIMPORANG115 + region_out
        os.makedirs(folder_out, exist_ok=True)
        
        if args.substitutionmodel == False:
            label = "singlent_filtered"
        else:
            label = "trinucleotide_filtered"
        
        file_out = folder_out + "/" + label + ".fa"
        noncds_alns.write(file_out)
        
        with open(folder_out + "/" + label + "_alnstat.txt", mode = "w") as out: 
            out.write("alignment length: " + str(len(noncds_alns)))

if __name__ == "__main__":
    main()

