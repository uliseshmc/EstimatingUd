"""Alignment-quality branch (exploratory): repeat-masked and sp_score-filtered intergenicAR.

Five steps, each one CLI command under ``ud quality``:

masked_databank  <region>/alldata_chrm<chr> -> <region>/nonmasked_alldata_chrm<chr>/ (one file per
                 fragment, renamed, degenerate columns dropped with gaps allowed, columns gapped
                 in all three species removed; fragments that fail go to the store's not_completed)
masked_concat    the same filter, concatenated into <region>/chrm<chr>/singlent_nomasked.fa
filter_spscore   nonmasked_alldata_chrm<chr> -> quality_aligned_chrm<chr>/
                 trinucleotide_filtered_window<W>_qualitythreshold<T>.fa
fit_spscore      GT_CpG_ss fit of that file -> quality_aligned_chrm<chr>/sm_output/
                 trinucleotide_lh_window<W>_qualitythreshold<T>.pickle
ens_spscore      Human-edge ENS of every fit_spscore pickle, chr 1-22 + X ->
                 results/tables/trinuc_ENS_spscore_window<W>_qualitythreshold<T>.csv
"""

from __future__ import annotations

import pickle
from pathlib import Path

import cogent3
import pandas as pd
from cogent3 import get_app, load_aligned_seqs
from cogent3.maths.matrix_exponential_integration import expected_number_subs

from estimating_ud import apps
from estimating_ud.regions import (
    CHROMOSOMES_NO_Y,
    MASKED_REGIONS,
    NEUTRAL_REGION,
    QUALITY_REGIONS,
    QUALITY_THRESHOLD,
    QUALITY_WINDOW,
    raw_dir,
    region_dir,
)
from estimating_ud.stages._common import concat_completed, open_fasta_store
from estimating_ud.stages.constraint import human_edge
from estimating_ud.stages.fit import trinuc_model_noncds

GAP_ALLOWED_FRAC = 0.90
ENS_COLUMNS = ["Region", "Chromosome", "window", "threshold", "ENS", "aln_length"]


def nonmasked_dir(data_root: Path, region: str, chrm: str) -> Path:
    return Path(data_root) / region / f"nonmasked_alldata_chrm{chrm}"


def quality_dir(data_root: Path, region: str, chrm: str) -> Path:
    return Path(data_root) / region / f"quality_aligned_chrm{chrm}"


def quality_fasta(data_root: Path, region: str, chrm: str, window: int, threshold: float) -> Path:
    return quality_dir(data_root, region, chrm) / f"trinucleotide_filtered_window{window}_qualitythreshold{threshold}.fa"


def quality_pickle(data_root: Path, region: str, chrm: str, window: int, threshold: float) -> Path:
    return quality_dir(data_root, region, chrm) / "sm_output" / f"trinucleotide_lh_window{window}_qualitythreshold{threshold}.pickle"


def ens_csv(tables_dir: Path, window: int, threshold: float) -> Path:
    return Path(tables_dir) / f"trinuc_ENS_spscore_window{window}_qualitythreshold{threshold}.csv"


def masked_app():
    loader = get_app("load_aligned", moltype="dna")
    rename = apps.renamer_noncds_aligned()
    omit_degs = get_app("omit_degenerates", gap_is_degen=False, moltype="dna", motif_length=1)
    remove_gaps = get_app("omit_gap_pos", allowed_frac=GAP_ALLOWED_FRAC)
    return loader + rename + omit_degs + remove_gaps


def masked_databank(chrm: str, data_root: Path, regions=tuple(MASKED_REGIONS)) -> list[Path]:
    written = []
    for region in regions:
        in_dstore = open_fasta_store(raw_dir(data_root, region, chrm))
        folder_out = nonmasked_dir(data_root, region, chrm)
        folder_out.mkdir(parents=True, exist_ok=True)
        out_dstore = cogent3.open_data_store(str(folder_out), suffix="fa", mode="w")
        app = masked_app() + get_app("write_seqs", data_store=out_dstore)
        _ = list(app.apply_to(in_dstore[:], parallel=False))
        written.append(folder_out)
    return written


def masked_concat(chrm: str, data_root: Path, regions=tuple(MASKED_REGIONS)) -> list[Path]:
    app = masked_app()
    written = []
    for region in regions:
        alns = concat_completed(app, open_fasta_store(raw_dir(data_root, region, chrm)))
        folder_out = region_dir(data_root, region, chrm)
        folder_out.mkdir(parents=True, exist_ok=True)
        file_out = folder_out / "singlent_nomasked.fa"
        alns.write(str(file_out))
        (folder_out / "singlent_nomasked_alnstat.txt").write_text(f"alignment length: {len(alns)}")
        written.append(file_out)
    return written


def spscore_app(window: int, threshold: float):
    # inputs are already renamed (Human/Chimpanzee/Orangutan) by masked_databank
    loader = get_app("load_aligned", moltype="dna")
    quality = apps.filter_windows_spscore(window=window, threshold_filter=threshold)
    omit_degs = get_app("omit_degenerates", moltype="dna", motif_length=3)
    return loader + quality + omit_degs


def filter_spscore(
    chrm: str, data_root: Path, window: int = QUALITY_WINDOW, threshold: float = QUALITY_THRESHOLD, regions=tuple(QUALITY_REGIONS)
) -> list[Path]:
    written = []
    for region in regions:
        alns = concat_completed(spscore_app(window, threshold), open_fasta_store(nonmasked_dir(data_root, region, chrm)))
        file_out = quality_fasta(data_root, region, chrm, window, threshold)
        file_out.parent.mkdir(parents=True, exist_ok=True)
        alns.write(str(file_out))
        written.append(file_out)
    return written


def fit_spscore(
    chrm: str, data_root: Path, window: int = QUALITY_WINDOW, threshold: float = QUALITY_THRESHOLD, regions=tuple(QUALITY_REGIONS)
) -> list[Path]:
    written = []
    for region in regions:
        aln = load_aligned_seqs(str(quality_fasta(data_root, region, chrm, window, threshold)), moltype="dna")
        result = trinuc_model_noncds()(aln)
        file_out = quality_pickle(data_root, region, chrm, window, threshold)
        file_out.parent.mkdir(parents=True, exist_ok=True)
        with open(file_out, "wb") as out:
            out.write(pickle.dumps(result))
        written.append(file_out)
    return written


def ens_spscore(
    data_root: Path,
    window: int = QUALITY_WINDOW,
    threshold: float = QUALITY_THRESHOLD,
    chromosomes=tuple(CHROMOSOMES_NO_Y),
    regions=tuple(QUALITY_REGIONS),
) -> pd.DataFrame:
    """Human-edge ENS of every ``fit_spscore`` pickle, one row per region and chromosome.

    Same layout as the stage 5 ENS table (intergenicAR labelled ``IGAR``, ``Chromosome`` as
    text) plus the ``window`` and ``threshold`` that select the pickle. ``aln_length`` is the
    length of the alignment kept inside the pickled result: no ``_alnstat.txt`` is written for
    the quality fasta. A missing pickle raises ``FileNotFoundError``; Y is excluded by default
    because no fragment of it keeps a window.
    """
    rows = []
    for region in regions:
        label = "IGAR" if region == NEUTRAL_REGION else region
        for chrm in chromosomes:
            with open(quality_pickle(data_root, region, chrm, window, threshold), "rb") as infile:
                result = pickle.load(infile)
            mprobs, Q = human_edge(result)
            rows.append(
                {
                    "Region": label,
                    "Chromosome": chrm,
                    "window": window,
                    "threshold": threshold,
                    "ENS": expected_number_subs(mprobs, Q, t=1.0),
                    "aln_length": len(result.alignment),
                }
            )
    return pd.DataFrame(rows, columns=ENS_COLUMNS)


def write_ens_spscore(
    data_root: Path,
    tables_dir: Path,
    window: int = QUALITY_WINDOW,
    threshold: float = QUALITY_THRESHOLD,
    chromosomes=tuple(CHROMOSOMES_NO_Y),
) -> Path:
    tables_dir = Path(tables_dir)
    tables_dir.mkdir(parents=True, exist_ok=True)
    path = ens_csv(tables_dir, window, threshold)
    ens_spscore(data_root, window, threshold, chromosomes).to_csv(path, index=False)
    return path
