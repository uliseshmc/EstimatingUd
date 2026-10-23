"""Helpers shared by the stages."""

from __future__ import annotations

from pathlib import Path

import cogent3
from cogent3 import get_app

from estimating_ud.regions import MOTIF_LENGTH, check_submodel


def omit_degenerates(submodel: str):
    """omit_degenerates with the motif length of the submodel (1 or 3)."""
    return get_app("omit_degenerates", moltype="dna", motif_length=MOTIF_LENGTH[check_submodel(submodel)])


def open_fasta_store(folder: Path):
    return cogent3.open_data_store(str(folder), suffix="fa", mode="r")


def concat_completed(app, in_dstore):
    """Apply ``app`` to every member of ``in_dstore`` and concatenate the completed results."""
    concat = get_app("concat", moltype="dna")
    completed = [r for r in app.as_completed(in_dstore[:], parallel=False) if r]
    return concat(completed)
