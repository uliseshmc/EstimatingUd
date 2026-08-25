from cogent3 import get_app
from cogent3 import load_aligned_seqs
import argparse

from EstimatingUd import paths

from EstimatingUd import trinuc_models as trinucs

#Same as testing_trinuc_model but with optimise_motif_probs set to True, which is the default. This is to test that the motif probabilities are being optimized correctly.
def noncds_model(upperparms, max_restarts, tolerance):
    sm_noncds=trinucs.GT_CpG_ss()
    paramnames = sm_noncds.get_param_list()
    rules_noncds = [{"par_name": n, "is_independent": True, "upper": upperparms} for n in paramnames if n!="omega"] + [{"par_name": "omega", "value": 1.0, "is_constant": True}]
    GT_subsmodel = get_app("model", "GT_CpG_ss",
                      show_progress=True,
                      optimise_motif_probs=True,
                      param_rules=rules_noncds,
                      opt_args={"max_restarts": max_restarts, "tolerance": tolerance}
                      )
    return GT_subsmodel

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-chrm", "--chromosome", type=str, required=True, help='Chromosome to process')
    parser.add_argument("-upper", "--upperparm", type=float, required=True, help='Upper bound for parameters in non-CDS regions')
    parser.add_argument("-maxr", "--max_restarts", type=int, required=True, help='Maximum number of restarts for optimization')
    parser.add_argument("-tol", "--tolerance", type=float, required=True, help='Tolerance for optimization')
    args = parser.parse_args()
    
    file_in = paths.DATA_HUMCHIMPORANG115 + "intergenicAR/chrm" + args.chromosome + "/trinucleotide_filtered.fa"
    aln = load_aligned_seqs(file_in, moltype="dna")

    GT_subsmodel = noncds_model(args.upperparm, args.max_restarts, args.tolerance)

    result = GT_subsmodel(aln)

    file_out = ("lfresults/motiffit_intergenicAR_chrm" + args.chromosome 
                + "upper" + str(args.upperparm) 
                + "_maxr" + str(args.max_restarts) 
                + "_tol" + str(args.tolerance) + ".txt")

    with open(file_out, "w") as f:
        print(result.lf, file=f)

if __name__ == "__main__":
    main()
