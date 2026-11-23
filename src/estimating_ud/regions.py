"""Region, chromosome and submodel vocabularies plus the on-disk layout under the data root.

The data layout is unchanged from the script generation of the pipeline:

    <data_root>/<region>/alldata_chrm<chr>/            raw eti output, one .fa per gene / block
    <data_root>/<region>/chrm<chr>/<label>.fa          concatenated filtered alignment
    <data_root>/<region>/chrm<chr>/<label>_alnstat.txt "alignment length: N"
    <data_root>/<region>/chrm<chr>/human_seq.fa        concatenated Human sequence
    <data_root>/<region>/chrm<chr>/sm_output/<fit>_lh.pickle

The three intron regions all read from ``introns/alldata_chrm<chr>`` and differ only in
the sampler applied.
"""

from __future__ import annotations

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
    "introns3UTR",
    "introns5UTR",
    "introns_nonUTR",
    "intergenicAR",
    "distalIG",
    "proximal5IG",
    "proximal3IG",
]
INTRON_REGIONS = ["introns3UTR", "introns5UTR", "introns_nonUTR"]
NONCDS_REGIONS = ["intergenicAR", "intronsAR", "distalIG", "proximal5IG", "proximal3IG"]
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
    """Directory of raw eti output for a region; intron regions share ``introns/``."""
    folder = "introns" if region in INTRON_REGIONS else region
    return Path(data_root) / folder / f"alldata_chrm{chrm}"


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
