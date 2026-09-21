#!/usr/bin/env python3
import cogent3
from cogent3 import get_app
import paths
import libs
import argparse
import os

CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y"]
REGIONS = ["cds", "intronsAR", "introns3UTR", "introns5UTR", "introns_nonUTR", "intergenicAR", "distalIG", "proximal5IG", "proximal3IG"]

def gethumanseqs_cds(mut_motif_length):
    loader = get_app("load_unaligned", moltype="dna")
    get_human_seq = libs.gethumanseq_cds_unaligned(motif_length=mut_motif_length)
    
    cds_app = loader + get_human_seq

    return cds_app

def gethumanseqs_noncds(mut_motif_length):
    loader = get_app("load_aligned", moltype="dna")
    get_human_seq = libs.gethumanseq_noncds_aligned(motif_length=mut_motif_length)
    
    noncds_app = loader + get_human_seq

    return noncds_app

def gethumanseqs_introns(region, mut_motif_length):
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

    get_human_seq = libs.gethumanseq_noncds_aligned(motif_length=mut_motif_length)

    introns_app = loader + rename_noncds + get_region + get_human_seq
    return introns_app

def get_humanseq_app(genomic_region):
    if genomic_region == "cds":
        app = gethumanseqs_cds(3)
    elif (genomic_region == "introns3UTR") | (genomic_region == "introns5UTR") | (genomic_region == "introns_nonUTR"):
        app = gethumanseqs_introns(genomic_region, 3)
    else:
        app = gethumanseqs_noncds(3)

    return app

def get_folder_in(genomic_region, chrom):
    if (genomic_region == "introns3UTR") | (genomic_region == "introns5UTR") | (genomic_region == "introns_nonUTR"):
        folder_in = "introns/alldata_chrm" + chrom
    else:
        folder_in = genomic_region + "/alldata_chrm" + chrom

    return folder_in
    

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-mutmotif",
        "--mutationmotiflength",
        type=int,
        required=True,
        help="Length of the motifs for the mutation model",
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

    for genomic_region in REGIONS:
        app = get_humanseq_app(genomic_region)
        
        relative_folder_in = get_folder_in(genomic_region, args.chromosome)
        folder_in = paths.DATA_HUMCHIMPORANGOR114 + relative_folder_in
        in_dstore = cogent3.open_data_store(folder_in, suffix='fa', mode='r')
        
        nonconcat_alns = [r.obj for r in app.as_completed(in_dstore[:], parallel=False) if r]
        concat_alns = cogent3.make_seq("".join(str(s) for s in nonconcat_alns), name="Human", moltype="dna")

        relative_folder_out = genomic_region + "/chrm" + args.chromosome
        folder_out = paths.DATA_HUMCHIMPORANGOR114 + relative_folder_out
        os.makedirs(folder_out, exist_ok=True)
        
        file_out = folder_out + "/human_seq_motiflength" + str(args.mutationmotiflength) + ".fa"
        concat_alns.write(file_out)
        
        with open(folder_out + "/human_seq_motiflength" + str(args.mutationmotiflength) + "_alnstat.txt", mode = "w") as out: 
            out.write("alignment length: " + str(len(concat_alns)))

if __name__ == "__main__":
    main()

