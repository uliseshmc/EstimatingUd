import cogent3
from cogent3 import get_app
from cogent3 import load_aligned_seqs
import paths
import pickle
import os


REGIONS = ["cds", "introns", "intergenicAR", "intronsAR", "distalIG", "proximal5IG", "proximal3IG"]
CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y"]

def fitting_cds():
    #Setting up the rules for model fitting of cds regions
    GN_subsmodel = get_app("model", "GN", time_het="max", lf_args={"discrete_edges": ["Gorilla"]}, optimise_motif_probs=True, show_progress=True)

    return GN_subsmodel

def fitting_noncds():
    GN_subsmodel = get_app("model", "GN", time_het="max", lf_args={"discrete_edges": ["Gorilla"]}, optimise_motif_probs=True, show_progress=True)

    return GN_subsmodel

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
    region = args.region + "/alldata_chrm" + args.chromosome
    folder_in = paths.DATA_HUMCHIMPORANG115 + region
    file_in = folder_in + "filtered.fa"
    alns = load_aligned_seqs(file_in, moltype="dna")
    alns.source = args.region + "_chrm" + args.chromosome

    data_out = folder_in + "submodels_output/"
    os.makedirs(data_out, exist_ok=True)
    with open(data_out + "singlent_aln_stats.txt", "w") as f:
        f.write("alignment length: " + str(len(alns)))

    if args.region == "cds":
        subsmodel = fitting_cds()

    else:
        subsmodel = fitting_noncds()
    
    result_sm = subsmodel(alns)
    with open(data_out + "singlent_likelihood.pickle", mode = "wb") as out: 
        out.write(pickle.dumps(result_sm))


if __name__ == "__main__":
    main()