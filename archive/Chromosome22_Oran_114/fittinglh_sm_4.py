import argparse
import cogent3
from cogent3 import get_app
from cogent3 import load_aligned_seqs
import paths
import pickle
import os
import trinuc_models as trinucs # this module must be in the same directory as this notebook

SUBMODELS = ["trinuc", "singlent"]
REGIONS = ["cds", "intronsAR", "introns3UTR", "introns5UTR", "introns_nonUTR", "intergenicAR", "distalIG", "proximal5IG", "proximal3IG"]

def singlentmodel_cds():
    GN_subsmodel = get_app("model", "GN", time_het="max", lf_args={"discrete_edges": ["Orangutan"]}, optimise_motif_probs=False, show_progress=False)
    return GN_subsmodel

def singlentmodel_noncds():
    GN_subsmodel = get_app("model", "GN", time_het="max", lf_args={"discrete_edges": ["Orangutan"]}, optimise_motif_probs=False, show_progress=False)
    return GN_subsmodel

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

def submodel_cds(substitutionmodel):
    if substitutionmodel == "singlent":
        submodel = singlentmodel_cds()
    elif substitutionmodel == "trinuc":
        submodel = trinucmodel_cds()
    else:
        raise ValueError("Trying to use a substitution model other than singlent, or trinuc")
    return submodel

def submodel_noncds(substitutionmodel):
    if substitutionmodel == "singlent":
        submodel = singlentmodel_noncds()
    elif substitutionmodel == "trinuc":
        submodel = trinucmodel_noncds()
    else:
        raise ValueError("Trying to use a substitution model other than singlent, or trinuc")
    return submodel

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

    for genomic_region in REGIONS:

        relative_folder_in = genomic_region + "/chrm22" 
        folder_in = paths.DATA_HUMCHIMPORANGOR114 + relative_folder_in

        if args.substitutionmodel == "singlent":
            file_in = folder_in + "/singlent_filtered.fa"
            label_out = "/singlent"
        elif args.substitutionmodel == "trinuc":
            file_in = folder_in + "/trinucleotide_filtered.fa"
            label_out = "/trinucleotide"
        else:
            raise ValueError("Trying to use a substitution model other than singlent, or trinuc")
        
        alns = load_aligned_seqs(file_in, moltype="dna")

        if genomic_region == "cds":
            sm = submodel_cds(args.substitutionmodel)
        else:
            sm = submodel_noncds(args.substitutionmodel)

        result_sm = sm(alns)

        folder_out = folder_in + "/sm_output"
        os.makedirs(folder_out, exist_ok=True)

        with open(folder_out + label_out + "_lh.pickle", mode = "wb") as out: 
            out.write(pickle.dumps(result_sm))

if __name__ == "__main__":
    main()
