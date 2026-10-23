"""Stage 5: expected number of substitutions (ENS) on the Human edge and constraint.

For every region and chromosome the ENS is computed under the region's own fitted rate
matrix and under the ``intergenicAR`` matrix of the same chromosome (``selfQ`` 1 and 0 in
the ENS table). ``constraint = (ENS_IGAR - ENS) / ENS_IGAR``. ``intergenicAR`` itself appears
in the ENS table as ``IGAR``.
"""

from __future__ import annotations

import pickle
import re
from pathlib import Path

import cogent3
import numpy as np
import pandas as pd
from cogent3.maths.matrix_exponential_integration import expected_number_subs

from estimating_ud.regions import (
    CHROMOSOMES_NO_Y,
    CONSTRAINT_REGIONS,
    NEUTRAL_REGION,
    RESULTS_LABEL,
    SUBMODELS,
    alnstat_path,
    fit_pickle,
)

STOP_CODON_PSEUDOCOUNT = 1e-12


def parse_alignment_length(s: str, fmt: str = "{:.6e}") -> str:
    """Parse "alignment length: X" and return X formatted (scientific notation by default)."""
    m = re.search(r"alignment\s*length\s*:\s*([0-9]+(?:\.[0-9]+)?)", s, re.I)
    if not m:
        msg = f"Can't parse alignment length from: {s!r}"
        raise ValueError(msg)
    return fmt.format(float(m.group(1)))


def get_aln_length(data_root: Path, region: str, chrm: str, submodel: str) -> str:
    return parse_alignment_length(alnstat_path(data_root, region, chrm, submodel).read_text().strip())


def trinuc_alphabet():
    """All 64 trinucleotides, in cogent3 order; used to pad the stop codons into cds motif probs."""
    return cogent3.get_moltype("dna").alphabet.get_kmer_alphabet(3)


def load_fit(data_root: Path, region: str, chrm: str, submodel: str):
    with open(fit_pickle(data_root, region, chrm, submodel), "rb") as infile:
        return pickle.load(infile)


def human_edge(result):
    """(motif probabilities, uncalibrated Human-edge rate matrix) of a fitted model result."""
    mprobs = result.lf.get_param_value("mprobs")
    Q = result.lf.get_rate_matrix_for_edge("Human", calibrated=False)
    return mprobs, Q


def region_ens(result, region: str, submodel: str, Q_igar, alphabet):
    """ENS under the region's own Q and under the intergenicAR Q."""
    mprobs, Q = human_edge(result)
    ENS = expected_number_subs(mprobs, Q, t=1.0)

    if submodel == "trinuc" and region == "cds":
        # Trinucleotide cds models have no stop-codon states, but the IGAR matrix does:
        # pad the three stop codons in with a negligible probability.
        probs = result.lf.get_motif_probs().to_dict()
        freqs = np.array([probs.get(m, STOP_CODON_PSEUDOCOUNT) / (1 + 3 * STOP_CODON_PSEUDOCOUNT) for m in alphabet])
        ENS_igar = expected_number_subs(freqs, Q_igar, t=1.0)
    else:
        ENS_igar = expected_number_subs(mprobs, Q_igar, t=1.0)
    return ENS, ENS_igar


def calculate_constraint(ENS, ENS_neutral):
    return (ENS_neutral - ENS) / ENS_neutral


def compute(data_root: Path, submodel: str, chromosomes=tuple(CHROMOSOMES_NO_Y), regions=tuple(CONSTRAINT_REGIONS)):
    """Return (ENS table, constraint table) as DataFrames for one submodel."""
    alphabet = trinuc_alphabet()
    rows_ens, rows_constraint = [], []
    for chrm in chromosomes:
        igar = load_fit(data_root, NEUTRAL_REGION, chrm, submodel)
        mprobs_igar, Q_igar = human_edge(igar)
        rows_ens.append(
            {
                "Region": "IGAR",
                "Chromosome": chrm,
                "selfQ": 1,
                "ENS": expected_number_subs(mprobs_igar, Q_igar, t=1.0),
                "aln_length": get_aln_length(data_root, NEUTRAL_REGION, chrm, submodel),
            }
        )
        for region in regions:
            aln_length = get_aln_length(data_root, region, chrm, submodel)
            ENS, ENS_igar = region_ens(load_fit(data_root, region, chrm, submodel), region, submodel, Q_igar, alphabet)
            rows_ens.append({"Region": region, "Chromosome": chrm, "selfQ": 1, "ENS": ENS, "aln_length": aln_length})
            rows_ens.append({"Region": region, "Chromosome": chrm, "selfQ": 0, "ENS": ENS_igar, "aln_length": aln_length})
            rows_constraint.append({"Region": region, "Chromosome": chrm, "Constraint": calculate_constraint(ENS, ENS_igar)})

    ens = pd.DataFrame(rows_ens, columns=["Region", "Chromosome", "selfQ", "ENS", "aln_length"])
    constraint = pd.DataFrame(rows_constraint, columns=["Region", "Chromosome", "Constraint"])
    return ens, constraint


def run(data_root: Path, tables_dir: Path, submodels=tuple(SUBMODELS), chromosomes=tuple(CHROMOSOMES_NO_Y)) -> list[Path]:
    tables_dir = Path(tables_dir)
    tables_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for submodel in submodels:
        ens, constraint = compute(data_root, submodel, chromosomes)
        label = RESULTS_LABEL[submodel]
        for name, table in ((f"{label}_ENS.csv", ens), (f"{label}_constraint.csv", constraint)):
            path = tables_dir / name
            table.to_csv(path, index=False)
            written.append(path)
    return written
