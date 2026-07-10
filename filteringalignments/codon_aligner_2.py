#!/usr/bin/env python3
import cogent3
from cogent3 import get_app
from EstimatingUd import paths
from EstimatingUd import libs
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-chrm", "--chromosome", type=str, required=True, help='Chromosome to process')
    args = parser.parse_args()

    # set folders to read and export data
    folder_in = paths.DATA_HUMCHIMPORANG115 + 'cds/chrm' + args.chromosome + '/'
    folder_out = paths.DATA_HUMCHIMPORANG115 + 'cds/chrm' + args.chromosome + '/codon_aligned'

    in_dstore = cogent3.open_data_store(folder_in, suffix='fa', mode='r')
    out_dstore = cogent3.open_data_store(folder_out, suffix='fa', mode='w')

    print('Input data store:', in_dstore.describe)

    # perform a codon alignment and remove stop codons
    loader = get_app('load_unaligned', moltype='dna')
    rename = libs.renamer_cds_unaligned()
    trim_stops = get_app('trim_stop_codons')
    codon_align = get_app('progressive_align', 'codon', guide_tree='(Human:0.06,Chimpanzee:0.06,Orangutan:0.14)')
    writer = get_app('write_seqs', data_store=out_dstore)
    app = loader + rename + trim_stops + codon_align + writer

    # underscore is used to ignore the returned list object
    # By using this renamer function sequences with paralogs are removed
    _ = list(app.apply_to(in_dstore[:], parallel=False))

    print('Output data store:', out_dstore.describe)
    print('Output summary completed status:', out_dstore.summary_not_completed)


if __name__ == '__main__':
    main()
