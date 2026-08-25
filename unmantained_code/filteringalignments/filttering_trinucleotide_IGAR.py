#!/usr/bin/env python3
import cogent3
from cogent3 import get_app
from EstimatingUd import paths
from EstimatingUd import libs

def filtering_noncds():
    loader = get_app("load_aligned", moltype="dna")
    rename_noncds = libs.renamer_noncds_aligned()
    omit_gap_pos_app = get_app("omit_gap_pos", moltype="dna")
    omit_degs_noncds = get_app("omit_degenerates", moltype="dna", motif_length=3)
    noncds_app = loader + rename_noncds + omit_gap_pos_app + omit_degs_noncds

    return noncds_app

def main():
    noncds_app = filtering_noncds()
    concat = get_app("concat", moltype="dna")

    for i in range(1, 22):
        region = "intergenicAR/chrm" + str(i)
        folder_in = paths.DATA_HUMCHIMPORANG115 + region
        in_dstore = cogent3.open_data_store(folder_in, suffix='fa', mode='r')
        
        nonconcat_noncds = [r for r in noncds_app.as_completed(in_dstore[:], parallel=False) if r]
        noncds_alns = concat(nonconcat_noncds)

        file_out = paths.DATA_HUMCHIMPORANG115 + region + "/trinucleotide_filtered.fa"
        noncds_alns.write(file_out)

if __name__ == "__main__":
    main()
