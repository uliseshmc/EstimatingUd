import cogent3
import pandas as pd
import paths
import pickle
import re

import trinuc_models as trinucs # this module must be in the same directory as this notebook


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
    folder_in = paths.DATA_HUMCHIMPORANG115 + region_path

    if substitutionmodel == "singlent":
        file_in = folder_in + "trinucleotide_filtered_alnstat.txt"
    elif substitutionmodel == "trinuc":
        file_in = folder_in + "singlent_filtered_alnstat.txt"
    else:
        raise ValueError("Trying to use a substitution model other than singlent, or trinuc")
    
    with open(file_in, "r") as f:
        content = f.read().strip()
    alignment_length = parse_alignment_length(content)
    
    return alignment_length

def get_trinuc_alphabet():
    aln = cogent3.get_dataset("primate-brca1")
    aln = aln.take_seqs(["Human", "Chimpanzee", "Rhesus"])
    aln = aln.omit_gap_pos(allowed_gap_frac=0, motif_length=3)
    alphabet_trinucs = aln.moltype.alphabet.get_kmer_alphabet(3)
    
    return alphabet_trinucs

def get_model_path(region, chrm):
    region_path = region + "/chrm" + chrm + "/sm_output/"
    folder_in = paths.DATA_HUMCHIMPORANG115 + region_path

    if substitutionmodel == "singlent":
        file_in = folder_in + "trinucleotide_filtered_alnstat.txt"
    elif substitutionmodel == "trinuc":
        file_in = folder_in + "singlent_filtered_alnstat.txt"
    else:
        raise ValueError("Trying to use a substitution model other than singlent, or trinuc")
    
    file_in = folder_in + "trinuc_lh.pickle"

    return file_in

# Given a region and chromosome, return the ENS given a trinucleotide likelihood result
def get_region_ENS(region, chrm):

    file_in = get_model_path(region, chrm)
    with open(file_in, mode = "rb") as infile: 
        result_sm=pickle.load(infile)

    ENS_cpg = result_sm.lf.get_scaled_lengths("CpG")['Human']
    ENS_notcpg = result_sm.lf.get_scaled_lengths("notCpG")['Human']
    
    return ENS_cpg, ENS_notcpg

# Given a region and chromosome, return the ENS calculated from replacing Q with the Q from IGAR alignments
def get_region_ENS_QIGAR(region, chrm, result_IGAR, alphabet_trinucs):
    if (region == "cds"):
        ENS_cpg_QIGAR, ENS_notcpg_QIGAR = get_region_ENS_QIGAR_cds(chrm, result_IGAR, alphabet_trinucs)
    else:
        ENS_cpg_QIGAR, ENS_notcpg_QIGAR = get_region_ENS_QIGAR_noncds(region, chrm, result_IGAR)
    return ENS_cpg_QIGAR, ENS_notcpg_QIGAR

def get_region_ENS_QIGAR_cds(chrm, result_IGAR, alphabet_trinucs):

    file_in = get_model_path("cds", chrm)
    with open(file_in, mode = "rb") as infile: 
        result_sm=pickle.load(infile)

    mprobs = result_sm.lf.get_motif_probs().to_dict()
    mle_freqs = {alphabet_trinuc: mprobs.get(alphabet_trinuc, 1e-12) / (1+3e-12) for alphabet_trinuc in alphabet_trinucs}

    new_lf = trinucs.modified_lf(result_IGAR.lf)
    new_lf.set_motif_probs(mle_freqs)

    ENS_cpg_QIGAR = new_lf.get_scaled_lengths("CpG")["Human"]
    ENS_notcpg_QIGAR = new_lf.get_scaled_lengths("notCpG")["Human"]
    
    return ENS_cpg_QIGAR, ENS_notcpg_QIGAR

def get_region_ENS_QIGAR_noncds(region, chrm, result_IGAR):
    file_in = get_model_path(region, chrm)
    with open(file_in, mode = "rb") as infile: 
        result_sm=pickle.load(infile)

    mprobs = result_sm.lf.get_motif_probs().to_dict()

    new_lf = trinucs.modified_lf(result_IGAR.lf)
    new_lf.set_motif_probs(mprobs)

    ENS_cpg_QIGAR = new_lf.get_scaled_lengths("CpG")["Human"]
    ENS_notcpg_QIGAR = new_lf.get_scaled_lengths("notCpG")["Human"]
    
    return ENS_cpg_QIGAR, ENS_notcpg_QIGAR

def calculate_constraint(ENS, ENSneutral):
    return (ENSneutral - ENS)/ENSneutral 

def main():
    # this alphabet is used in cds regions to ensure that the motif probabilities are calculated for all trinucleotides, 
    # including stop codons
    alphabet_trinucs = get_trinuc_alphabet()

    row_data_ENS = []
    row_data_constraint = []

    for chromosome in CHROMOSOMES:
        file_in = get_model_path("intergenicAR", chromosome)
        with open(file_in, mode = "rb") as infile: 
            result_IGAR=pickle.load(infile)

        aln_length = get_aln_length("intergenicAR", chromosome)
        IGAR_ENS_cpg = result_IGAR.lf.get_scaled_lengths("CpG")['Human']
        IGAR_ENS_notcpg = result_IGAR.lf.get_scaled_lengths("notCpG")['Human']

        row_data_ENS.append({
            "Region": "IGAR",
            "Chromosome": chromosome,
            "CpG or nonCpG": "CpG",
            "selfQ": 1,
            "ENS": IGAR_ENS_cpg,
            "aln length": aln_length
        })
        row_data_ENS.append({
            "Region": "IGAR",
            "Chromosome": chromosome,
            "CpG or nonCpG": "nonCpG",
            "selfQ": 1,
            "ENS": IGAR_ENS_notcpg,
            "aln length": aln_length
        })

        for region in REGIONS:

            aln_length = get_aln_length(region, chromosome)
            ENS_cpg, ENS_notcpg = get_region_ENS(region, chromosome)
            row_data_ENS.append({
                "Region": region,
                "Chromosome": chromosome,
                "CpG or nonCpG": "CpG",
                "selfQ": 1,
                "ENS": ENS_cpg,
                "aln length": aln_length
            })
            row_data_ENS.append({
                "Region": region,
                "Chromosome": chromosome,
                "CpG or nonCpG": "nonCpG",
                "selfQ": 1,
                "ENS": ENS_notcpg,
                "aln length": aln_length
            })

            ENS_cpg_QIGAR, ENS_notcpg_QIGAR = get_region_ENS_QIGAR(region, chromosome, result_IGAR, alphabet_trinucs)
            row_data_ENS.append({
                "Region": region,
                "Chromosome": chromosome,
                "CpG or nonCpG": "CpG",
                "selfQ": 0,
                "ENS": ENS_cpg_QIGAR,
                "aln length": aln_length
            })
            row_data_ENS.append({
                "Region": region,
                "Chromosome": chromosome,
                "CpG or nonCpG": "nonCpG",
                "selfQ": 0,
                "ENS": ENS_notcpg_QIGAR,
                "aln length": aln_length
            })

            cpg_constraint = calculate_constraint(ENS_cpg, ENS_cpg_QIGAR)
            notcpg_constraint = calculate_constraint(ENS_notcpg, ENS_notcpg_QIGAR)
            row_data_constraint.append({
                "Region": region,
                "Chromosome": chromosome,
                "CpG or nonCpG": "CpG",
                "Constraint": cpg_constraint
            })
            row_data_constraint.append({
                "Region": region,
                "Chromosome": chromosome,
                "CpG or nonCpG": "nonCpG",
                "Constraint": notcpg_constraint
            })

    #selfQ=1 for all regions, but 0 when using the Q from IGAR alignments
    ENS_data = pd.DataFrame(row_data_ENS, columns=["Region", "Chromosome", "CpG or nonCpG", "selfQ", "ENS", "aln length"])

    file_out = "output_data/trinuc_ENS.csv"
    ENS_data.to_csv(file_out, index=False)

    constraint_data = pd.DataFrame(row_data_constraint, columns=["Region", "Chromosome", "CpG or nonCpG", "Constraint"])
    file_out = "output_data/trinuc_constraint.csv"
    constraint_data.to_csv(file_out, index=False)

    stats = constraint_data.groupby(['Region', 'CpG or nonCpG'])['Constraint'].agg(['mean', 'sem']).unstack()
    file_out = "output_data/trinuc_constraint_mean_sem.csv"
    stats.to_csv(file_out, index=False)

if __name__ == "__main__":
    main()