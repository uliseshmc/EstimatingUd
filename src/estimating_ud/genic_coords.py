"""Where the 5' UTR, 3' UTR and intron classes are, from the installed gene annotation.

The classes are defined by coordinates, not by where the mask character falls in an
alignment (docs/todo.md, "The UTR split does not delimit UTRs"). Decisions taken 2026-11-23:

1. UTR = exons of a gene's canonical transcript outside that transcript's CDS; 5' or 3' by
   the gene's strand.
2. One intron class = protein-coding gene spans minus canonical-transcript exons, with
   repeats masked.
3. Where annotations of different genes overlap, a position goes to cds first, then UTR,
   then intron, so every position is counted once.

Every interval is 0-based and half-open, ``[start, stop)``: the convention of the
``gene_attr`` / ``transcript_attr`` tables written by ``eti install`` and of the
``--ref_coords`` files read by ``eti alignments``.

Repeats are NOT removed here. These tables say where each class lies; repeat positions inside
them are masked when the alignments are extracted (``eti alignments --mask <repeat list>``),
as for the proximal and distal intergenic regions. Sizes computed from these tables therefore
include repeat positions and are upper bounds on what the extracted, repeat-masked classes hold.

Positions that are 5' UTR of one gene and 3' UTR of another are kept apart as
``utr_ambiguous`` and belong to neither UTR class (decision of 2026-11-23: leave them out;
0.023 Mb genome-wide, 0.06% of all UTR).

Write the per-chromosome tables with::

    ud genic-coords              # all chromosomes
    ud genic-coords --chrm 22    # one chromosome

(``python -m estimating_ud.genic_coords`` does the same), which writes ``<data root>/genic_coordinates/chrom<chr>_{utr5,utr3,genes_merged}.tsv``.
"""

from __future__ import annotations

import dataclasses
import io
from pathlib import Path
from typing import Iterable, Sequence

Interval = tuple[int, int]

CLASSES = ("cds", "utr5", "utr3", "utr_ambiguous", "exon_other", "introns")
COORD_TABLES = ("utr5", "utr3", "genes_merged")  # the tables eti alignments needs
COORD_DIR = "genic_coordinates"
REF_COORDS_HEADER = ("species", "seqid", "start", "stop", "strand")

GENE_SQL = """
SELECT ga.stable_id, ga.symbol, ga.seqid, ga.start, ga.stop, ga.strand,
       ta.transcript_stable_id, ta.transcript_biotype, ta.transcript_spans, ta.cds_spans
FROM gene_attr ga
LEFT JOIN transcript_attr ta ON ta.transcript_id = ga.canonical_transcript_id
WHERE ga.biotype = ? AND ga.seqid = ?
ORDER BY ga.start, ga.stop, ga.stable_id
"""


# ---------------------------------------------------------------- interval arithmetic


def merge(intervals: Iterable[Interval]) -> list[Interval]:
    """Sorted, non-overlapping, non-touching union of intervals; empty intervals dropped."""
    out: list[Interval] = []
    for start, stop in sorted((int(a), int(b)) for a, b in intervals if b > a):
        if out and start <= out[-1][1]:
            if stop > out[-1][1]:
                out[-1] = (out[-1][0], stop)
        else:
            out.append((start, stop))
    return out


def subtract(a: Sequence[Interval], b: Sequence[Interval]) -> list[Interval]:
    """Parts of ``a`` not covered by ``b``. Both must already be merged."""
    out: list[Interval] = []
    j = 0
    for start, stop in a:
        cur = start
        while j < len(b) and b[j][1] <= cur:
            j += 1
        k = j
        while k < len(b) and b[k][0] < stop:
            if b[k][0] > cur:
                out.append((cur, b[k][0]))
            cur = max(cur, b[k][1])
            k += 1
        if cur < stop:
            out.append((cur, stop))
    return out


def intersect(a: Sequence[Interval], b: Sequence[Interval]) -> list[Interval]:
    """Parts covered by both ``a`` and ``b``. Both must already be merged."""
    out: list[Interval] = []
    i = j = 0
    while i < len(a) and j < len(b):
        start, stop = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if start < stop:
            out.append((start, stop))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def total(intervals: Iterable[Interval]) -> int:
    return sum(stop - start for start, stop in intervals)


# ---------------------------------------------------------------- one gene


@dataclasses.dataclass(frozen=True)
class Gene:
    """A protein-coding gene with the exons and CDS of its canonical transcript."""

    stable_id: str
    seqid: str
    start: int
    stop: int
    strand: int
    exons: tuple[Interval, ...] = ()
    cds: tuple[Interval, ...] = ()
    symbol: str = ""
    transcript_id: str = ""
    transcript_biotype: str = ""


@dataclasses.dataclass(frozen=True)
class GeneParts:
    utr5: tuple[Interval, ...]
    utr3: tuple[Interval, ...]
    other: tuple[Interval, ...]  # exon sequence that is neither CDS nor UTR (see gene_parts)
    cds_outside_exons: int  # should be 0; counted so a broken annotation row is visible


def gene_parts(gene: Gene) -> GeneParts:
    """Split one gene's canonical exons into 5' UTR, 3' UTR and the rest.

    Exon sequence below the lowest CDS position is 5' UTR on the plus strand and 3' UTR on
    the minus strand; above the highest CDS position it is the reverse. A transcript with
    no CDS has no UTR: all its exons go to ``other``. So does non-CDS exon sequence that
    lies between the first and last CDS position.
    """
    if gene.strand not in (1, -1):
        msg = f"{gene.stable_id}: strand must be 1 or -1, got {gene.strand!r}"
        raise ValueError(msg)
    exons, cds = merge(gene.exons), merge(gene.cds)
    if not cds:
        return GeneParts((), (), tuple(exons), 0)
    low_edge, high_edge = cds[0][0], cds[-1][1]
    low, high, other = [], [], []
    for start, stop in subtract(exons, cds):
        if stop <= low_edge:
            low.append((start, stop))
        elif start >= high_edge:
            high.append((start, stop))
        else:
            other.append((start, stop))
    utr5, utr3 = (low, high) if gene.strand == 1 else (high, low)
    return GeneParts(tuple(utr5), tuple(utr3), tuple(other), total(subtract(cds, exons)))


# ---------------------------------------------------------------- one chromosome


@dataclasses.dataclass
class Partition:
    """Disjoint classes covering the merged protein-coding gene spans of one chromosome."""

    cds: list[Interval]
    utr5: list[Interval]
    utr3: list[Interval]
    utr_ambiguous: list[Interval]
    exon_other: list[Interval]
    introns: list[Interval]
    genes_merged: list[Interval]
    utr5_plus: list[Interval]  # raw 5' UTR of plus-strand genes, for the strand column
    utr5_minus: list[Interval]
    utr3_plus: list[Interval]
    utr3_minus: list[Interval]
    raw: dict[str, int]  # sizes before the precedence rule, summed over genes (overlaps counted twice)

    def classes(self) -> dict[str, list[Interval]]:
        return {name: getattr(self, name) for name in CLASSES}

    def sizes(self) -> dict[str, int]:
        return {name: total(ivs) for name, ivs in self.classes().items()}


def build_partition(genes: Iterable[Gene]) -> Partition:
    """Apply cds > UTR > intron across all genes of one chromosome."""
    genes = list(genes)
    seqids = {g.seqid for g in genes}
    if len(seqids) > 1:
        msg = f"genes from several chromosomes in one partition: {sorted(seqids)}"
        raise ValueError(msg)
    parts = [gene_parts(g) for g in genes]
    cds = merge(iv for g in genes for iv in g.cds)
    exons = merge(iv for g in genes for iv in g.exons)
    spans = merge((min(g.start, g.stop), max(g.start, g.stop)) for g in genes)

    def side(attr: str, strand: int) -> list[Interval]:
        return merge(iv for g, p in zip(genes, parts) if g.strand == strand for iv in getattr(p, attr))

    utr5_plus, utr5_minus = side("utr5", 1), side("utr5", -1)
    utr3_plus, utr3_minus = side("utr3", 1), side("utr3", -1)
    utr5 = subtract(merge(utr5_plus + utr5_minus), cds)
    utr3 = subtract(merge(utr3_plus + utr3_minus), cds)
    ambiguous = intersect(utr5, utr3)
    utr5, utr3 = subtract(utr5, ambiguous), subtract(utr3, ambiguous)
    assigned = merge(cds + utr5 + utr3 + ambiguous)
    # exons can reach past the recorded gene span in a broken row; keep the classes inside the spans
    exon_other = intersect(subtract(exons, assigned), spans)
    introns = subtract(spans, merge(assigned + exon_other))
    cds, utr5, utr3, ambiguous = (intersect(x, spans) for x in (cds, utr5, utr3, ambiguous))
    raw = {
        "cds": sum(total(merge(g.cds)) for g in genes),
        "utr5": sum(total(p.utr5) for p in parts),
        "utr3": sum(total(p.utr3) for p in parts),
        "exon_other": sum(total(p.other) for p in parts),
        "gene_span": sum(abs(g.stop - g.start) for g in genes),
        "cds_outside_exons": sum(p.cds_outside_exons for p in parts),
    }
    return Partition(cds, utr5, utr3, ambiguous, exon_other, introns, spans,
                     utr5_plus, utr5_minus, utr3_plus, utr3_minus, raw)


def overlaps_between(classes: dict[str, Sequence[Interval]]) -> dict[tuple[str, str], int]:
    """Positions shared by each pair of classes, found by one sweep over all interval ends.

    Deliberately does not use merge/subtract/intersect, so it can check their output.
    """
    events = []
    for name, intervals in classes.items():
        for start, stop in intervals:
            events.append((start, 1, name))
            events.append((stop, 0, name))
    events.sort()  # at equal position, closings (0) come before openings (1)
    shared: dict[tuple[str, str], int] = {}
    open_names: dict[str, int] = {}
    previous = None
    for position, opening, name in events:
        if previous is not None and position > previous:
            active = sorted(n for n, count in open_names.items() if count > 0)
            for i, first in enumerate(active):
                for second in active[i + 1:]:
                    shared[(first, second)] = shared.get((first, second), 0) + position - previous
                if open_names[first] > 1:
                    shared[(first, first)] = shared.get((first, first), 0) + position - previous
        open_names[name] = open_names.get(name, 0) + (1 if opening else -1)
        previous = position
    return shared


def strand_of(interval: Interval, plus: Sequence[Interval], minus: Sequence[Interval]) -> int:
    """1 or -1 when genes of one strand contribute the interval, 0 when both do."""
    on_plus = bool(intersect([interval], plus))
    on_minus = bool(intersect([interval], minus))
    if on_plus and on_minus:
        return 0
    return -1 if on_minus else 1


def ref_coords_rows(part: Partition, name: str, seqid: str, species: str = "homo_sapiens") -> list[tuple]:
    """Rows for an ``eti alignments --ref_coords`` file: species, seqid, start, stop, strand."""
    intervals = part.genes_merged if name == "genes_merged" else getattr(part, name)
    if name == "utr5":
        strands = [strand_of(iv, part.utr5_plus, part.utr5_minus) for iv in intervals]
    elif name == "utr3":
        strands = [strand_of(iv, part.utr3_plus, part.utr3_minus) for iv in intervals]
    else:
        strands = [1] * len(intervals)  # eti does not use the column; other regions also write 1
    return [(species, seqid, start, stop, strand) for (start, stop), strand in zip(intervals, strands)]


def write_ref_coords(path: Path, rows: Iterable[tuple], header: bool = True) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="\n") as out:
        if header:
            out.write("\t".join(REF_COORDS_HEADER) + "\n")
        for row in rows:
            out.write("\t".join(str(value) for value in row) + "\n")


# ---------------------------------------------------------------- reading the install


def _spans(blob) -> tuple[Interval, ...]:
    if blob is None:
        return ()
    import numpy

    array = blob if isinstance(blob, numpy.ndarray) else numpy.load(io.BytesIO(bytes(blob)))
    return tuple((int(start), int(stop)) for start, stop in array.reshape(-1, 2))


def genes_from_connection(conn, seqid: str, biotype: str = "protein_coding") -> list[Gene]:
    """Genes of one chromosome from a duckdb connection holding gene_attr and transcript_attr."""
    genes = []
    for (stable_id, symbol, chrom, start, stop, strand,
         transcript_id, transcript_biotype, exon_blob, cds_blob) in conn.sql(GENE_SQL, params=(biotype, str(seqid))).fetchall():
        genes.append(Gene(
            stable_id=str(stable_id), seqid=str(chrom),
            start=min(int(start), int(stop)), stop=max(int(start), int(stop)),
            strand=-1 if int(strand) < 0 else 1,
            exons=_spans(exon_blob), cds=_spans(cds_blob),
            symbol=symbol or "", transcript_id=transcript_id or "", transcript_biotype=transcript_biotype or "",
        ))
    return genes


def open_annotation(install: Path, species: str = "human"):
    """(duckdb connection, genome name, installed config) for one species of an eti install."""
    from ensembl_tui import _config as eti_config
    from ensembl_tui import _genome as eti_genome

    config = eti_config.read_installed_cfg(install)
    genome_name = config.species_map.get_genome_name(species)
    annot = eti_genome.load_annotations_for_species(path=config.installed_genome(species=genome_name))
    return annot.genes.conn, genome_name, config


# ---------------------------------------------------------------- writing the tables


def coord_table_path(outdir: Path, seqid: str, name: str) -> Path:
    return Path(outdir) / f"chrom{seqid}_{name}.tsv"


def write_chromosome(conn, seqid: str, outdir: Path, species: str = "homo_sapiens") -> dict[str, int]:
    """Write the ``--ref_coords`` tables of one chromosome; return the size of every class.

    ``utr5`` and ``utr3`` are extracted as they are (with repeats masked). ``genes_merged`` holds
    the merged gene spans, from which the intron class is extracted by masking exons and repeats.
    """
    genes = genes_from_connection(conn, seqid)
    if not genes:
        msg = f"no protein-coding genes found on chromosome {seqid!r}"
        raise ValueError(msg)
    part = build_partition(genes)
    shared = overlaps_between(part.classes())
    if shared:
        msg = f"chromosome {seqid}: classes overlap, refusing to write tables: {shared}"
        raise RuntimeError(msg)
    for name in COORD_TABLES:
        write_ref_coords(coord_table_path(outdir, seqid, name), ref_coords_rows(part, name, str(seqid), species))
    return part.sizes()


def run(data_root: Path, chromosomes: Iterable[str], outdir: Path | None = None, report=None) -> list[Path]:
    """Write the coordinate tables of every chromosome in ``chromosomes``; return the paths written.

    ``report``, if given, is called with one line of class sizes per chromosome and a total line.
    """
    data_root = Path(data_root)
    outdir = Path(outdir) if outdir else data_root / COORD_DIR
    conn, genome_name, _ = open_annotation(data_root / "install")
    totals = dict.fromkeys(CLASSES, 0)
    written = []
    for seqid in chromosomes:
        sizes = write_chromosome(conn, seqid, outdir, species=genome_name)
        for name, size in sizes.items():
            totals[name] += size
        if report:
            report(f"chromosome {seqid}: " + ", ".join(f"{name} {size}" for name, size in sizes.items()))
        written.extend(coord_table_path(outdir, seqid, name) for name in COORD_TABLES)
    if report:
        report("total: " + ", ".join(f"{name} {size}" for name, size in totals.items()))
    return written


def main(argv: Sequence[str] | None = None) -> int:
    import argparse

    from estimating_ud.config import data_root
    from estimating_ud.regions import CHROMOSOMES

    parser = argparse.ArgumentParser(
        prog="python -m estimating_ud.genic_coords",
        description="Write the 5' UTR, 3' UTR and merged gene span coordinate tables read by eti alignments --ref_coords.",
    )
    parser.add_argument("--chrm", action="append", help="chromosome; repeat for several (default: 1-22, X, Y)")
    parser.add_argument("--data-root", default=None, help="alignment download root (default: DATA_HUMCHIMPORANGOR114)")
    parser.add_argument("--outdir", default=None, help=f"default: <data root>/{COORD_DIR}")
    args = parser.parse_args(argv)

    written = run(data_root(args.data_root), args.chrm or CHROMOSOMES, args.outdir, report=lambda line: print(line, flush=True))
    print(f"tables written to {written[0].parent}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
