"""Command-line interface: ``ud <stage> --chrm <chr> [--submodel singlent|trinuc]``.

Each command only parses options and calls the matching ``run`` function in
:mod:`estimating_ud.stages`. The stage modules are imported lazily so that ``ud --help``
does not pay for importing cogent3.
"""

from __future__ import annotations

from pathlib import Path

import click

from estimating_ud.config import DATA_ENV, RESULTS_ENV, results_dir
from estimating_ud.regions import CHROMOSOMES, CHROMOSOMES_NO_Y, QUALITY_THRESHOLD, QUALITY_WINDOW, SUBMODELS


def chrm_option(f):
    # "-chrm" is kept as an alias so the pre-restructure docs and wrappers stay valid
    return click.option(
        "-chrm", "--chrm", "chrm", required=True, type=click.Choice(CHROMOSOMES), help="Chromosome to process."
    )(f)


def submodel_option(f):
    return click.option(
        "-submodel",
        "--submodel",
        "submodel",
        required=True,
        type=click.Choice(SUBMODELS),
        help="Substitution model family: singlent (motif length 1) or trinuc (motif length 3).",
    )(f)


def data_root_option(f):
    return click.option(
        "--data-root",
        "data_root",
        envvar=DATA_ENV,
        required=True,
        type=click.Path(file_okay=False, exists=True, path_type=Path),
        help=f"Alignment download root. Defaults to ${DATA_ENV}.",
    )(f)


def window_options(f):
    f = click.option(
        "--window", default=QUALITY_WINDOW, show_default=True, type=int, help="Window size in alignment columns."
    )(f)
    return click.option(
        "--threshold",
        default=QUALITY_THRESHOLD,
        show_default=True,
        type=float,
        help="Fraction of one species' window that must be non-gap for the threshold alignment.",
    )(f)


@click.group()
@click.version_option(package_name="estimating-ud")
def main():
    """Estimate Ud from Human/Chimpanzee/Orangutan alignments, one stage at a time."""


@main.command("align-cds")
@chrm_option
@submodel_option
@data_root_option
def align_cds(chrm, submodel, data_root):
    """Stage 1: codon-align cds homologs, trim stops, drop degenerate motifs."""
    from estimating_ud.stages import align_cds as stage

    click.echo(stage.run(chrm, submodel, data_root))


@main.command("filter-introns")
@chrm_option
@submodel_option
@data_root_option
def filter_introns(chrm, submodel, data_root):
    """Stage 2: sample introns5UTR, introns3UTR and introns_nonUTR from the intron alignments."""
    from estimating_ud.stages import filter_introns as stage

    for path in stage.run(chrm, submodel, data_root):
        click.echo(path)


@main.command("filter-noncds")
@chrm_option
@submodel_option
@data_root_option
def filter_noncds(chrm, submodel, data_root):
    """Stage 3: filter intergenicAR, intronsAR, distalIG, proximal5IG, proximal3IG."""
    from estimating_ud.stages import filter_noncds as stage

    for path in stage.run(chrm, submodel, data_root):
        click.echo(path)


@main.command("human-seq")
@chrm_option
@data_root_option
def human_seq(chrm, data_root):
    """Concatenate the Human sequence of all nine regions (writes human_seq.fa)."""
    from estimating_ud.stages import human_seq as stage

    for path in stage.run(chrm, data_root):
        click.echo(path)


@main.command("fit")
@chrm_option
@submodel_option
@data_root_option
def fit(chrm, submodel, data_root):
    """Stage 4: fit the substitution model of every region and pickle the result."""
    from estimating_ud.stages import fit as stage

    for path in stage.run(chrm, submodel, data_root):
        click.echo(path)


@main.command("constraint")
@data_root_option
@click.option(
    "--submodel",
    "submodels",
    multiple=True,
    type=click.Choice(SUBMODELS),
    default=SUBMODELS,
    show_default=True,
    help="Submodel(s) to summarise; repeatable.",
)
@click.option(
    "--chrm",
    "chromosomes",
    multiple=True,
    type=click.Choice(CHROMOSOMES_NO_Y),
    default=CHROMOSOMES_NO_Y,
    help="Chromosome(s) to include; repeatable. Default: 1-22 and X.",
)
@click.option(
    "--results-dir",
    "results",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help=f"Where to write the tables. Defaults to <repo>/results/tables or ${RESULTS_ENV}/tables.",
)
def constraint(data_root, submodels, chromosomes, results):
    """Stage 5: ENS and constraint tables for every region (writes results/tables/*.csv)."""
    from estimating_ud.stages import constraint as stage

    tables_dir = (Path(results) / "tables") if results else results_dir("tables")
    for path in stage.run(data_root, tables_dir, submodels, chromosomes):
        click.echo(path)


@main.group()
def quality():
    """Alignment-quality branch: repeat-masked and sp_score-filtered intergenicAR."""


@quality.command("masked-databank")
@chrm_option
@data_root_option
def masked_databank(chrm, data_root):
    """Write one renamed, gap-filtered file per fragment to <region>/nonmasked_alldata_chrm<chr>."""
    from estimating_ud.stages import quality as stage

    for path in stage.masked_databank(chrm, data_root):
        click.echo(path)


@quality.command("masked-concat")
@chrm_option
@data_root_option
def masked_concat(chrm, data_root):
    """Concatenate the gap-filtered fragments into <region>/chrm<chr>/singlent_nomasked.fa."""
    from estimating_ud.stages import quality as stage

    for path in stage.masked_concat(chrm, data_root):
        click.echo(path)


@quality.command("filter-spscore")
@chrm_option
@data_root_option
@window_options
def filter_spscore(chrm, data_root, window, threshold):
    """Keep the windows whose sp_score beats the threshold alignment; concatenate them."""
    from estimating_ud.stages import quality as stage

    for path in stage.filter_spscore(chrm, data_root, window, threshold):
        click.echo(path)


@quality.command("fit-spscore")
@chrm_option
@data_root_option
@window_options
def fit_spscore(chrm, data_root, window, threshold):
    """Fit GT_CpG_ss to the sp_score-filtered alignment and pickle the result."""
    from estimating_ud.stages import quality as stage

    for path in stage.fit_spscore(chrm, data_root, window, threshold):
        click.echo(path)


@quality.command("ens-spscore")
@data_root_option
@window_options
@click.option(
    "--chrm",
    "chromosomes",
    multiple=True,
    type=click.Choice(CHROMOSOMES_NO_Y),
    default=CHROMOSOMES_NO_Y,
    help="Chromosome(s) to include; repeatable. Default: 1-22 and X.",
)
@click.option(
    "--results-dir",
    "results",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help=f"Where to write the table. Defaults to <repo>/results/tables or ${RESULTS_ENV}/tables.",
)
def ens_spscore(data_root, window, threshold, chromosomes, results):
    """ENS of every fit-spscore pickle (writes results/tables/trinuc_ENS_spscore_*.csv)."""
    from estimating_ud.stages import quality as stage

    tables_dir = (Path(results) / "tables") if results else results_dir("tables")
    click.echo(stage.write_ens_spscore(data_root, tables_dir, window, threshold, chromosomes))
