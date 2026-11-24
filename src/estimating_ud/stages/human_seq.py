"""Concatenate the Human sequence of every region (input to the mutability notebook)."""

from __future__ import annotations

from pathlib import Path

import cogent3
from cogent3 import get_app

from estimating_ud import apps
from estimating_ud.regions import REGIONS, check_coordinate_named, human_seq_fasta, raw_dir
from estimating_ud.stages._common import open_fasta_store


def human_seq_app(region: str):
    if region == "cds":
        return get_app("load_unaligned", moltype="dna") + apps.gethumanseq_cds_unaligned()
    return get_app("load_aligned", moltype="dna") + apps.gethumanseq_noncds_aligned()


def run(chrm: str, data_root: Path, regions=tuple(REGIONS)) -> list[Path]:
    written = []
    for region in regions:
        app = human_seq_app(region)
        check_coordinate_named(raw_dir(data_root, region, chrm), region)
        in_dstore = open_fasta_store(raw_dir(data_root, region, chrm))
        pieces = [r.obj for r in app.as_completed(in_dstore[:], parallel=False) if r]
        human = cogent3.make_seq("".join(str(s) for s in pieces), name="Human", moltype="dna")

        file_out = human_seq_fasta(data_root, region, chrm)
        file_out.parent.mkdir(parents=True, exist_ok=True)
        human.write(str(file_out))
        (file_out.parent / "human_seq_alnstat.txt").write_text(f"alignment length: {len(human)}")
        written.append(file_out)
    return written
