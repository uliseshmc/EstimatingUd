import argparse
import cogent3
from cogent3 import get_app
from cogent3 import load_aligned_seqs
import paths
import pickle
import os
import trinuc_models as trinucs # this module must be in the same directory as this notebook

#pipline only works for trinuc and intergenicAR
CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y"]
REGIONS = ["intergenicAR"]


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
        "-chrm",
        "--chromosome",
        type=str,
        required=True,
        choices=CHROMOSOMES,
        help="Chromosome to process (one of: " + ", ".join(CHROMOSOMES) + ")",
    )
    parser.add_argument("-window", 
        type=int, 
        default=60, 
        help="window size in alignment columns (default 60)")
    parser.add_argument(
        "-threshold",
        type=float,
        default=0.5,
        help="fraction of one species' window that must be non-gap for the threshold alignment (default 0.5)",
    )
    args = parser.parse_args()

    for genomic_region in REGIONS:

        relative_folder_in = genomic_region + "/quality_aligned_chrm" + args.chromosome + "/"
        folder_in = paths.DATA_HUMCHIMPORANGOR114 + relative_folder_in
        file_in = folder_in + "trinucleotide_filtered_window" + str(args.window) + "_qualitythreshold" + str(args.threshold) + ".fa"
        label_out = "/trinucleotide"
        
        alns = load_aligned_seqs(file_in, moltype="dna")
        sm = trinucmodel_noncds()
        result_sm = sm(alns)

        folder_out = folder_in + "/sm_output"
        os.makedirs(folder_out, exist_ok=True)

        with open(folder_out + label_out + "_lh.pickle", mode = "wb") as out: 
            out.write(pickle.dumps(result_sm))

if __name__ == "__main__":
    main()
