"""Stage 1: codon-align the cds homolog collections, trim stop codons, drop degenerate motifs."""

from __future__ import annotations

from pathlib import Path

import numpy
from cogent3 import get_app

from estimating_ud import apps
from estimating_ud.regions import alnstat_path, filtered_fasta, raw_dir, write_alnstat
from estimating_ud.stages._common import concat_completed, omit_degenerates, open_fasta_store

# cogent3.maths.util sets numpy.seterr(divide="raise") at import time, while
# cogent3.evolve.likelihood_tree sets seterr(all="ignore") -- whichever loads last
# wins. progressive_align's pair-HMM always takes log(0) on the structurally-zero
# BEGIN/END entries of its transition matrix (a legitimate -inf), so if "raise" is
# in effect every alignment fails with FloatingPointError. Pin it here.
numpy.seterr(divide="ignore")

GUIDE_TREE = "(Human:0.06,Chimpanzee:0.06,Orangutan:0.14)"


def cds_app(submodel: str):
    loader = get_app("load_unaligned", moltype="dna")
    rename = apps.renamer_cds_unaligned()
    trim_stops = get_app("trim_stop_codons")
    codon_align = get_app("progressive_align", "codon", guide_tree=GUIDE_TREE)
    return loader + rename + trim_stops + codon_align + omit_degenerates(submodel)


def run(chrm: str, submodel: str, data_root: Path) -> Path:
    in_dstore = open_fasta_store(raw_dir(data_root, "cds", chrm))
    alns = concat_completed(cds_app(submodel), in_dstore)

    file_out = filtered_fasta(data_root, "cds", chrm, submodel)
    file_out.parent.mkdir(parents=True, exist_ok=True)
    alns.write(str(file_out))
    write_alnstat(alnstat_path(data_root, "cds", chrm, submodel), len(alns))
    return file_out
