#This code was run using the environment cog2026.10.22_En0.7.9
import cogent3
from cogent3 import get_app
import paths
import libs
import argparse
import os

CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y"]
REGIONS = ["intergenicAR", "intronsAR", "distalIG"]

def filtering_noncds(out_dstore):
    loader = get_app("load_aligned", moltype="dna")
    rename_noncds = libs.renamer_noncds_aligned()
    omit_degs = get_app("omit_degenerates", gap_is_degen=False, moltype="dna", motif_length=1)
    remove_gaps = get_app("omit_gap_pos", allowed_frac=0.90)
    writer = get_app('write_seqs', data_store=out_dstore)

    noncds_app = loader + rename_noncds + omit_degs + remove_gaps + writer

    return noncds_app

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

    for genomic_region in REGIONS:
        relative_folder_in = genomic_region + "/alldata_chrm" + args.chromosome
        folder_in = paths.DATA_HUMCHIMPORANGOR114 + relative_folder_in
        in_dstore = cogent3.open_data_store(folder_in, suffix='fa', mode='r')

        relative_folder_out = genomic_region + "/nonmasked_alldata_chrm" + args.chromosome
        folder_out = paths.DATA_HUMCHIMPORANGOR114 + relative_folder_out

        os.makedirs(folder_out, exist_ok=True)
        out_dstore = cogent3.open_data_store(folder_out, suffix='fa', mode='w')

        noncds_app = filtering_noncds(out_dstore)
        _ = list(noncds_app.apply_to(in_dstore[:], parallel=False))

if __name__ == "__main__":
    main()

