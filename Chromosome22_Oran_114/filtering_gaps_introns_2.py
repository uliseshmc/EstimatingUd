#!/usr/bin/env python3
import cogent3
from cogent3 import get_app
import whole_genome_Orangutan.paths as paths
import whole_genome_Orangutan.libs as libs
import argparse
import os

SUBMODELS = ["trinuc", "singlent"]
INTRONREGIONS = ["introns3UTR", "introns5UTR", "introns_nonUTR"]

def filtering_introns(region, substitutionmodel):
    loader = get_app("load_aligned", moltype="dna")
    rename_noncds = libs.renamer_noncds_aligned()
    if region == "introns5UTR":
        get_region = libs.sample_UTR5()
    elif region == "introns3UTR":
        get_region = libs.sample_UTR3()
    elif region == "introns_nonUTR":
        get_region = libs.removeUTRs_fromintrons()
    else:
        raise ValueError("Trying to filter an intron region other than introns3UTR, introns5UTR, or introns_nonUTR")

    if substitutionmodel == "singlent":
        omit_degs = get_app("omit_degenerates", moltype="dna", motif_length=1)
    elif substitutionmodel == "trinuc":
        omit_degs = get_app("omit_degenerates", moltype="dna", motif_length=3)
    else:
        raise ValueError("Trying to use a substitution model other than singlent, or trinuc")

    introns_app = loader + rename_noncds + get_region + omit_degs
    return introns_app

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

    for region in INTRONREGIONS:
        
        introns_app = filtering_introns(region, args.substitutionmodel)
            
        region = "introns/alldata_chrm22"
        folder_in = paths.DATA_HUMCHIMPORANG115 + region
        in_dstore = cogent3.open_data_store(folder_in, suffix='fa', mode='r')
        
        nonconcat_introns = [r for r in introns_app.as_completed(in_dstore[:], parallel=False) if r]
        introns_alns = concat(nonconcat_introns)
    
        region_out = region + "/chrm22"
        folder_out = paths.DATA_HUMCHIMPORANG115 + region_out
        os.makedirs(folder_out, exist_ok=True)

        if args.substitutionmodel == "singlent":
            label = "singlent_filtered"
        elif substitutionmodel == "trinuc":
            label = "trinucleotide_filtered"
        else:
            raise ValueError("Trying to use a substitution model other than singlent, or trinuc")
            
        file_out = folder_out + "/" + label + ".fa"
        introns_alns.write(file_out)

        with open(folder_out + "/" + label + "_alnstat.txt", mode = "w") as out: 
            out.write("alignment length: " + str(len(introns_alns)))

if __name__ == "__main__":
    main()

