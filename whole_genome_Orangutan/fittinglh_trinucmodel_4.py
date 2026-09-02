import argparse
import cogent3
from cogent3 import get_app
from cogent3 import load_aligned_seqs
import matplotlib.pyplot as plt
import paths
import pickle
import os

import trinuc_models as trinucs # this module must be in the same directory as this notebook

REGIONS = ["cds", "introns", "introns3UTR", "introns5UTR", "introns_nonUTR", "intergenicAR", "intronsAR", "distalIG", "proximal5IG", "proximal3IG"]
CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y"]

def trinucmodel_cds():
    sm_noncds=trinucs.GNC_CpG_ss()
    paramnames = sm_noncds.get_param_list()
    rules_cds = [{"par_name": n, "is_independent": True} for n in paramnames]
    GNC_subsmodel = get_app("model", "GNC_CpG_ss",
                      show_progress=True,
                      param_rules=rules_cds)

    return GNC_subsmodel

def trinucmodel_noncds():
    #Setting up the rules for model fitting of noncds regions
    sm_noncds=trinucs.GT_CpG_ss()
    paramnames = sm_noncds.get_param_list()
    rules_noncds = [{"par_name": n, "is_independent": True} for n in paramnames if n!="omega"] + [{"par_name": "omega", "value": 1.0, "is_constant": True}]
    GT_subsmodel = get_app("model", "GT_CpG_ss",
                      show_progress=True,
                      optimise_motif_probs=False,
                      param_rules=rules_noncds)

    return GT_subsmodel

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

    region = args.region + "/chrm" + args.chromosome
    folder_in = paths.DATA_HUMCHIMPORANG115 + region
    file_in = folder_in + "/trinucleotide_filtered.fa"
    alns = load_aligned_seqs(file_in, moltype="dna")

    if args.region == "cds":
        sm = trinucmodel_cds()
        result_sm = sm(alns)
    else:
        sm = trinucmodel_noncds()
        result_sm = sm(alns)

    data_out = folder_in + "/sm_output"
    os.makedirs(data_out, exist_ok=True)

    with open(data_out + "/trinuc_lh.pickle", mode = "wb") as out: 
        out.write(pickle.dumps(result_sm))

    with open(data_out + "/trinuc_alnstat.txt", mode = "w") as out: 
        out.write("alignment length: " + str(len(alns)))

if __name__ == "__main__":
    main()
