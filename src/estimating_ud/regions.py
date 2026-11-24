"""Region, chromosome and submodel vocabularies plus the on-disk layout under the data root.

The data layout is unchanged from the script generation of the pipeline:

    <data_root>/<region>/alldata_chrm<chr>/            raw eti output, one .fa per gene / block
    <data_root>/<region>/chrm<chr>/<label>.fa          concatenated filtered alignment
    <data_root>/<region>/chrm<chr>/<label>_alnstat.txt "alignment length: N"
    <data_root>/<region>/chrm<chr>/human_seq.fa        concatenated Human sequence
    <data_root>/<region>/chrm<chr>/sm_output/<fit>_lh.pickle

``utr5``, ``utr3`` and ``introns`` are located by coordinates (``ud genic-coords``, see
:mod:`estimating_ud.genic_coords`) and extracted with ``eti alignments --ref_coords``, so their raw
fragments are named ``homo_sapiens-<chr>-<start>-<stop>[-<i>].fa``. Until 2026-11 the folder
``introns/alldata_chrm<chr>`` held one fragment per gene (``ENSG...fa``) from which three classes
were cut at the first and last masked position; :func:`check_coordinate_named` refuses that layout.
"""

from __future__ import annotations

import re
from pathlib import Path

CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y"]
CHROMOSOMES_NO_Y = CHROMOSOMES[:-1]  # constraint and mutability skip Y (no cds alignment)

SUBMODELS = ("singlent", "trinuc")
MOTIF_LENGTH = {"singlent": 1, "trinuc": 3}
# note: the trinuc label is "trinucleotide", not "trinuc", on disk
FILTERED_LABEL = {"singlent": "singlent_filtered", "trinuc": "trinucleotide_filtered"}
FIT_LABEL = {"singlent": "singlent", "trinuc": "trinucleotide"}
RESULTS_LABEL = {"singlent": "singlent", "trinuc": "trinuc"}

REGIONS = [
    "cds",
    "intronsAR",
    "utr3",
    "utr5",
    "introns",
    "intergenicAR",
    "distalIG",
    "proximal5IG",
    "proximal3IG",
]
# regions whose location comes from the genic_coordinates tables (ud genic-coords)
GENIC_REGIONS = ["utr5", "utr3", "introns"]
NONCDS_REGIONS = ["intergenicAR", "intronsAR", "distalIG", "proximal5IG", "proximal3IG", *GENIC_REGIONS]
# every region whose constraint is measured against intergenicAR
CONSTRAINT_REGIONS = [r for r in REGIONS if r != "intergenicAR"]
NEUTRAL_REGION = "intergenicAR"
MASKED_REGIONS = ["intergenicAR", "intronsAR", "distalIG"]
QUALITY_REGIONS = ["intergenicAR"]
# sp_score window (columns, keep it a multiple of 3) and threshold fraction; the q2/q3 SLURM
# wrappers hard-code the same values
QUALITY_WINDOW = 60
QUALITY_THRESHOLD = 0.5


def check_submodel(submodel: str) -> str:
    if submodel not in SUBMODELS:
        msg = f"Trying to use a substitution model other than {', '.join(SUBMODELS)}: {submodel!r}"
        raise ValueError(msg)
    return submodel


def raw_dir(data_root: Path, region: str, chrm: str) -> Path:
    """Directory of raw eti output for a region."""
    return Path(data_root) / region / f"alldata_chrm{chrm}"


_COORDINATE_NAMED = re.compile(r"^[A-Za-z_]+-[^-]+-\d+-\d+(-\d+)?$")


def check_coordinate_named(folder: Path, region: str) -> None:
    """Refuse a raw folder of a coordinate-defined region that holds fragments named any other way.

    ``eti alignments --ref_coords`` names its output ``<species>-<chr>-<start>-<stop>[-<i>].fa``.
    A file such as ``ENSG00000100345.fa`` in ``introns/alldata_chrm<chr>`` is the earlier per-gene
    download, where coding sequence and repeats were masked but UTR exons were not.
    """
    if region not in GENIC_REGIONS:
        return
    for path in sorted(Path(folder).glob("*.fa")):
        if not _COORDINATE_NAMED.match(path.stem):
            msg = (
                f"{folder} holds fragments that are not named by coordinates (for example {path.name}). "
                f"It looks like the per-gene download used before the {region!r} region was defined by "
                f"coordinates; extract it again with scripts/download/{region}.sh before filtering."
            )
            raise ValueError(msg)


def region_dir(data_root: Path, region: str, chrm: str) -> Path:
    return Path(data_root) / region / f"chrm{chrm}"


def filtered_fasta(data_root: Path, region: str, chrm: str, submodel: str) -> Path:
    return region_dir(data_root, region, chrm) / f"{FILTERED_LABEL[check_submodel(submodel)]}.fa"


def alnstat_path(data_root: Path, region: str, chrm: str, submodel: str) -> Path:
    return region_dir(data_root, region, chrm) / f"{FILTERED_LABEL[check_submodel(submodel)]}_alnstat.txt"


def human_seq_fasta(data_root: Path, region: str, chrm: str) -> Path:
    return region_dir(data_root, region, chrm) / "human_seq.fa"


def fit_pickle(data_root: Path, region: str, chrm: str, submodel: str) -> Path:
    return region_dir(data_root, region, chrm) / "sm_output" / f"{FIT_LABEL[check_submodel(submodel)]}_lh.pickle"


def write_alnstat(path: Path, length: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"alignment length: {length}")
