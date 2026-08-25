import cogent3
from cogent3 import get_app
from EstimatingUd import paths
from EstimatingUd import libs
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-chrm", "--chromosome", type=str, required=True, help='Chromosome to process')
args = parser.parse_args()

# Sampling 3' UTR
folder_in = paths.DATA_HUMCHIMPORANG115 + 'introns/chrm' + args.chromosome + '/'
folder_out = paths.DATA_HUMCHIMPORANG115 + 'introns/chrm' + args.chromosome + '/3UTR'

in_dstore = cogent3.open_data_store(folder_in, suffix='fa', mode='r')
out_dstore = cogent3.open_data_store(folder_out, suffix='fa', mode='w')

loader = get_app("load_aligned", moltype="dna")
rename = libs.renamer_noncds_aligned()
get_UTR3 = libs.sample_UTR3()
writer = get_app("write_seqs", data_store=out_dstore)
app_UTR3 = loader + rename + get_UTR3 + writer

# underscore is to specify a variable we are not gonna use later
# By using this renamer function I throw away sequences with paralogs
_ = list(app_UTR3.apply_to(in_dstore[:], parallel=False))

#Comment out this lines to see the output of the data store and the summary of not completed sequences
print(out_dstore.describe)
print(out_dstore.summary_not_completed)

# Sampling 5' UTR
folder_in = paths.DATA_HUMCHIMPORANG115 + 'introns/chrm' + args.chromosome + '/'
folder_out = paths.DATA_HUMCHIMPORANG115 + 'introns/chrm' + args.chromosome + '/5UTR'

in_dstore = cogent3.open_data_store(folder_in, suffix='fa', mode='r')
out_dstore = cogent3.open_data_store(folder_out, suffix='fa', mode='w')

loader = get_app("load_aligned", moltype="dna")
rename = libs.renamer_noncds_aligned()
get_UTR5 = libs.sample_UTR5()
writer = get_app("write_seqs", data_store=out_dstore)
app_UTR5 = loader + rename + get_UTR5 + writer

# underscore is to specify a variable we are not gonna use later
# By using this renamer function I throw away sequences with paralogs
_ = list(app_UTR5.apply_to(in_dstore[:], parallel=False))
print(out_dstore.describe)
out_dstore.summary_not_completed

# Removing UTRs from intron sequences
folder_in = paths.DATA_HUMCHIMPORANG115 + 'introns/chrm' + args.chromosome + '/'
folder_out = paths.DATA_HUMCHIMPORANG115 + 'introns/chrm' + args.chromosome + '/noUTRs'

in_dstore = cogent3.open_data_store(folder_in, suffix='fa', mode='r')
out_dstore = cogent3.open_data_store(folder_out, suffix='fa', mode='w')

loader = get_app("load_aligned", moltype="dna")
rename = libs.renamer_noncds_aligned()
remove_UTR = libs.removeUTRs_fromintrons()
writer = get_app("write_seqs", data_store=out_dstore)
app_noUTR = loader + rename + remove_UTR + writer

# underscore is to specify a variable we are not gonna use later
# By using this renamer function I throw away sequences with paralogs
_ = list(app_noUTR.apply_to(in_dstore[:], parallel=False))
print(out_dstore.describe)
print(out_dstore.summary_not_completed)
