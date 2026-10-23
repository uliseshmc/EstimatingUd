#!/usr/bin/env python3
# Run with the cog2026.10.22_En0.7.9 environment.
# Subsets the unmasked per-fragment alignments in <region>/nonmaksed_alldata_chrm<chr>
# (written by filtering_databank_masked_nointrons.py) to the windows whose sp_score is
# above the threshold alignment (libs.filter_windows_spscore) and writes one file per
# fragment, with the same name, to <region>/qualityaligned_alldata_chrm<chr>.
# Fragments with no passing window are recorded in the output store's not_completed.
import argparse
import os

import cogent3
from cogent3 import get_app

import libs
import paths

CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y"]
REGIONS = ["intergenicAR"]


def filtering_quality(window, threshold_filter):
    # inputs are already renamed (Human/Chimpanzee/Orangutan) and degenerate-filtered
    loader = get_app("load_aligned", moltype="dna")
    quality = libs.filter_windows_spscore(window=window, threshold_filter=threshold_filter)
    omit_degs = get_app("omit_degenerates", moltype="dna", motif_length=3)

    return loader + quality + omit_degs


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

    concat = get_app("concat", moltype="dna")

    for genomic_region in REGIONS:
        relative_folder_in = genomic_region + "/nonmasked_alldata_chrm" + args.chromosome
        folder_in = paths.DATA_HUMCHIMPORANGOR114 + relative_folder_in
        in_dstore = cogent3.open_data_store(folder_in, suffix="fa", mode="r")

        quality_app = filtering_quality(args.window, args.threshold)
        # in_dstore.completed skips the not_completed records left by filtering_databank_masked_nointrons.py
        nonconcat_alns = [r for r in quality_app.as_completed(in_dstore[:], parallel=False) if r]
        alns = concat(nonconcat_alns)

        relative_folder_out = genomic_region + "/quality_aligned_chrm" + args.chromosome + "/"
        folder_out = paths.DATA_HUMCHIMPORANGOR114 + relative_folder_out
        os.makedirs(folder_out, exist_ok=True)

        file_out = folder_out + "trinucleotide_filtered_window" + str(args.window) + "_qualitythreshold" + str(args.threshold) + ".fa"
        alns.write(file_out)       

if __name__ == "__main__":
    main()
