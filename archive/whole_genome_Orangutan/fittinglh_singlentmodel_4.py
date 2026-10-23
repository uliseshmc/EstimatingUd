import argparse
import cogent3
from cogent3 import get_app
from cogent3 import load_aligned_seqs
import matplotlib.pyplot as plt
import paths
import pickle
import os

REGIONS = ["cds", "introns", "introns3UTR", "introns5UTR", "introns_nonUTR", "intergenicAR", "intronsAR", "distalIG", "proximal5IG", "proximal3IG"]
CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y"]

def singlentmodel_cds():
    GN_subsmodel = get_app("model", "GN", time_het="max", lf_args={"discrete_edges": ["Orangutan"]}, optimise_motif_probs=False, show_progress=False)
    return GN_subsmodel

def singlentmodel_noncds():
    GN_subsmodel = get_app("model", "GN", time_het="max", lf_args={"discrete_edges": ["Orangutan"]}, optimise_motif_probs=False, show_progress=False)
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

    region = args.region + "/chrm" + args.chromosome
    folder_in = paths.DATA_HUMCHIMPORANG115 + region
    file_in = folder_in + "/filtered.fa"
    alns = load_aligned_seqs(file_in, moltype="dna")

    if args.region == "cds":
        sm = singlentmodel_cds()
        result_sm = sm(alns)
    else:
        sm = singlentmodel_noncds()
        result_sm = sm(alns)

    data_out = folder_in + "/sm_output"
    os.makedirs(data_out, exist_ok=True)

    with open(data_out + "/singlent_lh.pickle", mode = "wb") as out: 
        out.write(pickle.dumps(result_sm))

    with open(data_out + "/singlent_alnstat.txt", mode = "w") as out: 
        out.write("alignment length: " + str(len(alns)))

if __name__ == "__main__":
    main()
