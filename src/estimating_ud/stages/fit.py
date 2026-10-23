"""Stage 4: fit the substitution models and pickle the results.

Single-nucleotide: ``GN`` with ``time_het="max"`` and Orangutan as a discrete edge.
Trinucleotide cds: ``GNC_CpG_ss`` with every parameter independent.
Trinucleotide non-cds: ``GT_CpG_ss`` with ``omega`` fixed at 1.
"""

from __future__ import annotations

import pickle
from pathlib import Path

from cogent3 import get_app, load_aligned_seqs

from estimating_ud import models
from estimating_ud.regions import REGIONS, check_submodel, filtered_fasta, fit_pickle


def singlent_model():
    return get_app(
        "model",
        "GN",
        time_het="max",
        lf_args={"discrete_edges": ["Orangutan"]},
        optimise_motif_probs=False,
        show_progress=False,
    )


def trinuc_model_cds():
    paramnames = models.GNC_CpG_ss().get_param_list()
    rules = [{"par_name": n, "is_independent": True} for n in paramnames]
    return get_app("model", "GNC_CpG_ss", show_progress=True, param_rules=rules)


def trinuc_model_noncds():
    paramnames = models.GT_CpG_ss().get_param_list()
    rules = [{"par_name": n, "is_independent": True} for n in paramnames if n != "omega"]
    rules.append({"par_name": "omega", "value": 1.0, "is_constant": True})
    return get_app("model", "GT_CpG_ss", show_progress=True, optimise_motif_probs=False, param_rules=rules)


def model_for(region: str, submodel: str):
    check_submodel(submodel)
    if submodel == "singlent":
        return singlent_model()
    return trinuc_model_cds() if region == "cds" else trinuc_model_noncds()


def fit_alignment(aln, region: str, submodel: str):
    return model_for(region, submodel)(aln)


def run(chrm: str, submodel: str, data_root: Path, regions=tuple(REGIONS)) -> list[Path]:
    written = []
    for region in regions:
        aln = load_aligned_seqs(str(filtered_fasta(data_root, region, chrm, submodel)), moltype="dna")
        result = fit_alignment(aln, region, submodel)

        file_out = fit_pickle(data_root, region, chrm, submodel)
        file_out.parent.mkdir(parents=True, exist_ok=True)
        with open(file_out, "wb") as out:
            out.write(pickle.dumps(result))
        written.append(file_out)
    return written
