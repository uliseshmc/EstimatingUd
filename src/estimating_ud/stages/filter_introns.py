"""Stage 2: sample the 5' UTR, 3' UTR and non-UTR parts of the intron alignments."""

from __future__ import annotations

from pathlib import Path

from cogent3 import get_app

from estimating_ud import apps
from estimating_ud.regions import INTRON_REGIONS, alnstat_path, filtered_fasta, raw_dir, write_alnstat
from estimating_ud.stages._common import concat_completed, omit_degenerates, open_fasta_store

SAMPLERS = {
    "introns5UTR": apps.sample_UTR5,
    "introns3UTR": apps.sample_UTR3,
    "introns_nonUTR": apps.removeUTRs_fromintrons,
}


def introns_app(region: str, submodel: str):
    if region not in SAMPLERS:
        msg = f"Trying to filter an intron region other than {', '.join(SAMPLERS)}: {region!r}"
        raise ValueError(msg)
    loader = get_app("load_aligned", moltype="dna")
    rename = apps.renamer_noncds_aligned()
    return loader + rename + SAMPLERS[region]() + omit_degenerates(submodel)


def run(chrm: str, submodel: str, data_root: Path, regions=tuple(INTRON_REGIONS)) -> list[Path]:
    written = []
    for region in regions:
        in_dstore = open_fasta_store(raw_dir(data_root, region, chrm))
        alns = concat_completed(introns_app(region, submodel), in_dstore)

        file_out = filtered_fasta(data_root, region, chrm, submodel)
        file_out.parent.mkdir(parents=True, exist_ok=True)
        alns.write(str(file_out))
        write_alnstat(alnstat_path(data_root, region, chrm, submodel), len(alns))
        written.append(file_out)
    return written
