"""Stage 3: gap/degenerate filtering of the pre-aligned non-cds, non-intron regions."""

from __future__ import annotations

from pathlib import Path

from cogent3 import get_app

from estimating_ud import apps
from estimating_ud.regions import NONCDS_REGIONS, alnstat_path, filtered_fasta, raw_dir, write_alnstat
from estimating_ud.stages._common import concat_completed, omit_degenerates, open_fasta_store


def noncds_app(submodel: str):
    loader = get_app("load_aligned", moltype="dna")
    rename = apps.renamer_noncds_aligned()
    return loader + rename + omit_degenerates(submodel)


def run(chrm: str, submodel: str, data_root: Path, regions=tuple(NONCDS_REGIONS)) -> list[Path]:
    app = noncds_app(submodel)
    written = []
    for region in regions:
        in_dstore = open_fasta_store(raw_dir(data_root, region, chrm))
        alns = concat_completed(app, in_dstore)

        file_out = filtered_fasta(data_root, region, chrm, submodel)
        file_out.parent.mkdir(parents=True, exist_ok=True)
        alns.write(str(file_out))
        write_alnstat(alnstat_path(data_root, region, chrm, submodel), len(alns))
        written.append(file_out)
    return written
