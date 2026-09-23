import cogent3
from cogent3.maths.matrix_exponential_integration import expected_number_subs
import pandas as pd
import paths
import pickle
import re
import numpy as np

import trinuc_models as trinucs # this module must be in the same directory as this notebook

SUBMODELS = ["trinuc", "singlent"]
REGIONS = ["cds", "introns_nonUTR", "introns3UTR", "introns5UTR", "introns_nonUTR", "intronsAR", "distalIG", "proximal5IG", "proximal3IG"]
CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X"]

def parse_alignment_length(s: str, fmt: str = "{:.6e}") -> str:
    """
    Parse strings like "alignment length: X" and return X formatted in scientific notation.

    Parameters
    - s: input string containing a token like "alignment length: 12345"
    - fmt: a format string compatible with float.format, default "{:.6e}" (6 decimals in exp notation)

    Returns
    - str: the parsed number formatted according to `fmt` (scientific notation by default)

    Raises ValueError if the token is not found.
    """
    m = re.search(r'alignment\s*length\s*:\s*([0-9]+(?:\.[0-9]+)?)', s, re.I)
    if not m:
        raise ValueError(f"Can't parse alignment length from: {s!r}")
    token = m.group(1)
    # normalize to float for consistent formatting
    value = float(token)
    return fmt.format(value)

# Given a region and chromosome, return the alignment length in scientific notation
def get_aln_length(region, chrm, submodel):
    region_path = region + "/chrm" + chrm
    folder_in = paths.DATA_HUMCHIMPORANGOR114 + region_path

    if submodel == "singlent":
        file_in = folder_in + "/trinucleotide_filtered_alnstat.txt"
    elif submodel == "trinuc":
        file_in = folder_in + "/singlent_filtered_alnstat.txt"
    else:
        raise ValueError("Trying to use a substitution model other than singlent, or trinuc")
    
    with open(file_in, "r") as f:
        content = f.read().strip()
    alignment_length = parse_alignment_length(content)
    
    return alignment_length

def get_trinuc_alphabet():
    aln = cogent3.get_dataset("primate-brca1")
    alphabet_trinucs = aln.moltype.alphabet.get_kmer_alphabet(3)
    
    return alphabet_trinucs

def get_model_path(region, chrm, submodel):
    region_path = region + "/chrm" + chrm + "/sm_output"
    folder_in = paths.DATA_HUMCHIMPORANGOR114 + region_path

    if submodel == "singlent":
        file_in = folder_in + "/singlent_lh.pickle"
    elif submodel == "trinuc":
        file_in = folder_in + "/trinucleotide_lh.pickle"
    else:
        raise ValueError("Trying to use a substitution model other than singlent, or trinuc")

    return file_in

# Given a region and chromosome, return the ENS and the hypothetical ENS_QIGAR
def get_region_ENS_ENSQIGAR(region, chrm, submodel, QIGAR, alphabet_trinucs):

    file_in = get_model_path(region, chrm, submodel)
    with open(file_in, mode = "rb") as infile: 
        result_sm=pickle.load(infile)

    motif_probs = result_sm.lf.get_param_value("mprobs")
    Q = result_sm.lf.get_rate_matrix_for_edge("Human", calibrated=False)
    ENS = expected_number_subs(motif_probs, Q, t=1.0)

    if (submodel == "trinuc") & (region == "cds"):
        mprobs = result_sm.lf.get_motif_probs().to_dict()
        #Trinucleotide cds models ignore stop codons. 
        #Since they are present on the QIGAR matrix, I add them to the possible set of motifs with a negligible probability 1e-12.
        mle_freqs = {alphabet_trinuc: mprobs.get(alphabet_trinuc, 1e-12) / (1+3e-12) for alphabet_trinuc in alphabet_trinucs}
        mle_freqs = np.array(list(mle_freqs.values()))
        ENS_QIGAR = expected_number_subs(mle_freqs, QIGAR, t=1.0)
    else:
        ENS_QIGAR = expected_number_subs(motif_probs, QIGAR, t=1.0)
    
    return ENS, ENS_QIGAR

def calculate_constraint(ENS, ENSneutral):
    return (ENSneutral - ENS)/ENSneutral 

def main():
    # this alphabet is used in cds regions to ensure that the motif probabilities are calculated for all trinucleotides, 
    # including stop codons
    alphabet_trinucs = get_trinuc_alphabet()

    row_data_ENS = []
    row_data_constraint = []

    for submodel in SUBMODELS:
        for chromosome in CHROMOSOMES:
            file_in = get_model_path("intergenicAR", chromosome, submodel)
            with open(file_in, mode = "rb") as infile: 
                result_IGAR=pickle.load(infile)

            aln_length = get_aln_length("intergenicAR", chromosome, submodel)

            IGAR_motif_probs = result_IGAR.lf.get_param_value("mprobs")
            QIGAR = result_IGAR.lf.get_rate_matrix_for_edge("Human", calibrated=False)
            IGAR_ENS = expected_number_subs(IGAR_motif_probs, QIGAR, t=1.0)

            row_data_ENS.append({
                "Region": "IGAR",
                "Chromosome": chromosome,
                "selfQ": 1,
                "ENS": IGAR_ENS,
                "aln_length": aln_length
            })

            for region in REGIONS:

                aln_length = get_aln_length(region, chromosome, submodel)
                ENS, ENS_QIGAR = get_region_ENS_ENSQIGAR(region, chromosome, submodel, QIGAR, alphabet_trinucs)
                row_data_ENS.append({
                    "Region": region,
                    "Chromosome": chromosome,
                    "selfQ": 1,
                    "ENS": ENS,
                    "aln_length": aln_length
                })
                row_data_ENS.append({
                    "Region": region,
                    "Chromosome": chromosome,
                    "selfQ": 0,
                    "ENS": ENS_QIGAR,
                    "aln_length": aln_length
                })

                constraint_value = calculate_constraint(ENS, ENS_QIGAR)
                row_data_constraint.append({
                    "Region": region,
                    "Chromosome": chromosome,
                    "Constraint": constraint_value
                })

        #selfQ=1 for all regions, but 0 when using the Q from IGAR alignments
        ENS_data = pd.DataFrame(row_data_ENS, columns=["Region", "Chromosome", "selfQ", "ENS", "aln_length"])
        constraint_data = pd.DataFrame(row_data_constraint, columns=["Region", "Chromosome", "Constraint"])

        if submodel == "singlent":
            label_out = "singlent"
        elif submodel == "trinuc":
            label_out = "trinuc"
        else:
            raise ValueError("Trying to use a substitution model other than singlent, or trinuc")

        
        file_out = "output_data/" + label_out + "_ENS.csv"
        ENS_data.to_csv(file_out, index=False)

        file_out = "output_data/" + label_out + "_constraint.csv"
        constraint_data.to_csv(file_out, index=False)

        stats = constraint_data.groupby(['Region'])['Constraint'].agg(['mean', 'sem']).unstack()
        file_out = "output_data/" + label_out + "_constraint_mean_sem.csv"
        stats.to_csv(file_out, index=False)

if __name__ == "__main__":
    main()