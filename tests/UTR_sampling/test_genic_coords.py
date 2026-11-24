"""Tests of the coordinate-defined UTR and intron classes (synthetic genes, answers worked by hand)."""

import io
import random

import pytest

try:
    from estimating_ud import genic_coords as gc
except ImportError:  # standalone copy next to this file
    import genic_coords as gc

EXONS = ((100, 200), (300, 400), (500, 650))
CDS = ((150, 200), (300, 400), (500, 560))


def gene(name, strand, exons=EXONS, cds=CDS, start=None, stop=None, seqid="22"):
    lo = min(s for s, _ in exons) if start is None else start
    hi = max(e for _, e in exons) if stop is None else stop
    return gc.Gene(stable_id=name, seqid=seqid, start=lo, stop=hi, strand=strand, exons=tuple(exons), cds=tuple(cds))


def positions(intervals):
    return {p for start, stop in intervals for p in range(start, stop)}


# ---- interval arithmetic


def test_merge_joins_overlapping_and_touching_and_drops_empty():
    assert gc.merge([(5, 9), (1, 3), (3, 4), (8, 12), (20, 20)]) == [(1, 4), (5, 12)]


def test_subtract_and_intersect():
    a, b = [(0, 10), (20, 30)], [(2, 4), (8, 22), (29, 40)]
    assert gc.subtract(a, b) == [(0, 2), (4, 8), (22, 29)]
    assert gc.intersect(a, b) == [(2, 4), (8, 10), (20, 22), (29, 30)]
    assert gc.subtract(a, []) == a and gc.intersect(a, []) == []
    assert gc.total(a) == 20


def test_interval_arithmetic_against_position_sets():
    rng = random.Random(1)
    for _ in range(300):
        a = gc.merge((s, s + rng.randint(0, 9)) for s in (rng.randint(0, 80) for _ in range(rng.randint(0, 6))))
        b = gc.merge((s, s + rng.randint(0, 9)) for s in (rng.randint(0, 80) for _ in range(rng.randint(0, 6))))
        assert positions(gc.subtract(a, b)) == positions(a) - positions(b)
        assert positions(gc.intersect(a, b)) == positions(a) & positions(b)
        assert gc.merge(gc.subtract(a, b)) == gc.subtract(a, b)


# ---- one gene


def test_plus_strand_gene():
    parts = gc.gene_parts(gene("P", 1))
    assert parts.utr5 == ((100, 150),)
    assert parts.utr3 == ((560, 650),)
    assert parts.other == () and parts.cds_outside_exons == 0


def test_minus_strand_gene_swaps_the_ends():
    parts = gc.gene_parts(gene("M", -1))
    assert parts.utr5 == ((560, 650),)  # 5' UTR is at the high-coordinate end
    assert parts.utr3 == ((100, 150),)


def test_single_exon_gene():
    parts = gc.gene_parts(gene("S", 1, exons=((1000, 2000),), cds=((1200, 1700),)))
    assert parts.utr5 == ((1000, 1200),) and parts.utr3 == ((1700, 2000),)


def test_transcript_without_utr():
    parts = gc.gene_parts(gene("N", -1, exons=CDS, cds=CDS))
    assert parts.utr5 == () and parts.utr3 == ()


def test_utr_only_exons_and_exon_that_is_part_utr_part_cds():
    exons = ((0, 50), (100, 200), (300, 400), (450, 500))
    parts = gc.gene_parts(gene("U", 1, exons=exons, cds=((120, 200), (300, 380))))
    assert parts.utr5 == ((0, 50), (100, 120))
    assert parts.utr3 == ((380, 400), (450, 500))


def test_transcript_without_cds_has_no_utr():
    parts = gc.gene_parts(gene("X", 1, cds=()))
    assert parts.utr5 == () and parts.utr3 == () and parts.other == EXONS


def test_noncoding_exon_inside_the_coding_range_is_not_utr():
    parts = gc.gene_parts(gene("O", 1, exons=((0, 100), (200, 300), (400, 500)), cds=((50, 100), (400, 450))))
    assert parts.utr5 == ((0, 50),) and parts.utr3 == ((450, 500),) and parts.other == ((200, 300),)


def test_bad_strand_is_refused():
    with pytest.raises(ValueError):
        gc.gene_parts(gene("B", 0))


# ---- one chromosome


def test_single_gene_partition_sizes():
    part = gc.build_partition([gene("P", 1)])
    assert part.sizes() == {"cds": 210, "utr5": 50, "utr3": 90, "utr_ambiguous": 0, "exon_other": 0, "introns": 200}
    assert part.introns == [(200, 300), (400, 500)]
    assert sum(part.sizes().values()) == gc.total(part.genes_merged) == 550


def test_cds_of_another_gene_wins_over_utr():
    a = gene("A", 1)  # 3' UTR (560, 650)
    b = gene("B", -1, exons=((600, 700),), cds=((600, 700),))
    part = gc.build_partition([a, b])
    assert part.utr3 == [(560, 600)]
    assert (600, 700) in part.cds


def test_utr_of_one_gene_wins_over_intron_of_another():
    host = gene("H", 1, exons=((0, 100), (900, 1000)), cds=((0, 100), (900, 1000)))
    nested = gene("N", 1, exons=((400, 500),), cds=((420, 480),))
    part = gc.build_partition([host, nested])
    assert part.utr5 == [(400, 420)] and part.utr3 == [(480, 500)]
    assert part.introns == [(100, 400), (500, 900)]


def test_five_prime_and_three_prime_utr_overlap_is_kept_apart():
    a = gene("A", 1)  # 3' UTR (560, 650)
    c = gene("C", 1, exons=((620, 700), (800, 900)), cds=((680, 700), (800, 900)))  # 5' UTR (620, 680)
    part = gc.build_partition([a, c])
    assert part.utr_ambiguous == [(620, 650)]
    assert part.utr3 == [(560, 620)] and part.utr5 == [(100, 150), (650, 680)]


def test_gene_span_beyond_the_canonical_transcript_is_intron():
    part = gc.build_partition([gene("L", 1, start=0, stop=800)])
    assert part.introns == [(0, 100), (200, 300), (400, 500), (650, 800)]


def test_genes_of_two_chromosomes_are_refused():
    with pytest.raises(ValueError):
        gc.build_partition([gene("A", 1), gene("B", 1, seqid="21")])


def random_gene(rng, name):
    start = rng.randint(0, 600)
    n_exons = rng.randint(1, 5)
    cuts = sorted(rng.sample(range(start, start + 400), 2 * n_exons))
    exons = tuple((cuts[2 * i], cuts[2 * i + 1]) for i in range(n_exons))
    coding = sorted(positions(exons))
    if rng.random() < 0.15:
        cds = ()
    else:
        lo, hi = sorted(rng.sample(range(len(coding) + 1), 2))
        keep = set(coding[lo:hi])
        cds = tuple(gc.merge((p, p + 1) for p in keep)) if keep else ()
        cds = tuple(gc.intersect(list(cds), list(exons)))
    return gc.Gene(name, "22", start - rng.randint(0, 20), start + 400 + rng.randint(0, 20), rng.choice((1, -1)), exons, cds)


def test_partition_against_position_sets():
    """The classes are disjoint, cover the gene spans exactly, and follow cds > UTR > intron."""
    rng = random.Random(7)
    for trial in range(150):
        genes = [random_gene(rng, f"G{i}") for i in range(rng.randint(1, 6))]
        part = gc.build_partition(genes)
        span, cds, exon, u5, u3 = set(), set(), set(), set(), set()
        for g in genes:
            span |= positions([(g.start, g.stop)])
            exon |= positions(g.exons)
            cds |= positions(g.cds)
            if g.cds:
                lo, hi = min(positions(g.cds)), max(positions(g.cds))
                below = {p for p in positions(g.exons) - positions(g.cds) if p < lo}
                above = {p for p in positions(g.exons) - positions(g.cds) if p > hi}
                u5 |= below if g.strand == 1 else above
                u3 |= above if g.strand == 1 else below
        want = {
            "cds": cds,
            "utr5": u5 - cds - u3,
            "utr3": u3 - cds - u5,
            "utr_ambiguous": (u5 & u3) - cds,
            "exon_other": exon - cds - u5 - u3,
            "introns": span - exon,
        }
        got = {name: positions(ivs) for name, ivs in part.classes().items()}
        assert got == want, trial
        assert gc.overlaps_between(part.classes()) == {}
        assert sum(part.sizes().values()) == len(span) == gc.total(part.genes_merged)


def test_overlap_check_finds_overlaps():
    shared = gc.overlaps_between({"a": [(0, 10), (5, 8)], "b": [(8, 12)], "c": [(20, 30)]})
    assert shared == {("a", "a"): 3, ("a", "b"): 2}


# ---- files and the install


def test_ref_coords_rows_and_file(tmp_path):
    part = gc.build_partition([gene("P", 1), gene("M", -1, exons=((2000, 2100), (2200, 2300)), cds=((2050, 2100), (2200, 2250)))])
    rows = gc.ref_coords_rows(part, "utr5", "22")
    assert rows == [("homo_sapiens", "22", 100, 150, 1), ("homo_sapiens", "22", 2250, 2300, -1)]
    assert gc.ref_coords_rows(part, "genes_merged", "22") == [("homo_sapiens", "22", 100, 650, 1), ("homo_sapiens", "22", 2000, 2300, 1)]
    path = tmp_path / "coords" / "chrom22_utr5.tsv"
    gc.write_ref_coords(path, rows)
    lines = path.read_text().splitlines()
    assert lines[0] == "species\tseqid\tstart\tstop\tstrand"
    assert all(len(line.split("\t")) == 5 for line in lines)
    assert lines[2] == "homo_sapiens\t22\t2250\t2300\t-1"


def test_strand_column_is_zero_when_both_strands_contribute():
    plus = gene("P", 1, exons=((0, 100), (200, 300)), cds=((50, 100), (200, 300)))  # 5' UTR (0, 50)
    minus = gene("M", -1, exons=((10, 40),), cds=())  # no CDS: contributes no UTR
    both = gene("Q", -1, exons=((-200, -100), (20, 60)), cds=((-200, -100),))  # 5' UTR (20, 60)
    part = gc.build_partition([plus, minus, both])
    # (50, 60) is CDS of the plus gene, so it is not UTR; (0, 50) is 5' UTR of genes on both strands
    assert gc.ref_coords_rows(part, "utr5", "22") == [("homo_sapiens", "22", 0, 50, 0)]


def install_connection():
    """A duckdb connection with the two tables eti install writes, holding a few made-up genes."""
    duckdb = pytest.importorskip("duckdb")
    numpy = pytest.importorskip("numpy")

    def blob(spans):
        with io.BytesIO() as out:
            numpy.save(out, numpy.array(spans, dtype=numpy.int32))
            return out.getvalue()

    conn = duckdb.connect(":memory:")
    conn.sql("CREATE TABLE gene_attr (gene_id INTEGER, stable_id TEXT, biotype TEXT, canonical_transcript_id INTEGER, "
             "seqid TEXT, coord_system_name TEXT, start INTEGER, stop INTEGER, strand TINYINT, symbol TEXT, description TEXT)")
    conn.sql("CREATE TABLE transcript_attr (gene_id INTEGER, transcript_id INTEGER, seqid TEXT, coord_system_name TEXT, "
             "start INTEGER, stop INTEGER, strand TINYINT, transcript_spans BLOB, cds_spans BLOB, "
             "transcript_stable_id TEXT, transcript_biotype TEXT, cds_stable_id TEXT)")
    genes = [(1, "G1", "protein_coding", 11, "22", "chromosome", 100, 650, 1, "AAA", None),
             (2, "G2", "protein_coding", 21, "22", "chromosome", 90, 700, -1, None, None),
             (3, "G3", "lncRNA", 31, "22", "chromosome", 5, 50, 1, "LNC", None),
             (4, "G4", "protein_coding", 41, "21", "chromosome", 5, 50, 1, "OTHERCHR", None),
             (5, "G5", "protein_coding", 99, "22", "chromosome", 800, 900, 1, "NOTRANSCRIPT", None)]
    conn.executemany("INSERT INTO gene_attr VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", genes)
    transcripts = [(1, 11, "22", "chromosome", 100, 650, 1, blob(EXONS), blob(CDS), "T11", "protein_coding", "P11"),
                   (1, 12, "22", "chromosome", 100, 300, 1, blob([(100, 300)]), None, "T12", "retained_intron", None),
                   (2, 21, "22", "chromosome", 90, 700, -1, blob([(90, 700)]), None, "T21", "protein_coding", None),
                   (4, 41, "21", "chromosome", 5, 50, 1, blob([(5, 50)]), blob([(10, 40)]), "T41", "protein_coding", "P41")]
    conn.executemany("INSERT INTO transcript_attr VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", transcripts)
    return conn


def test_genes_from_connection_reads_the_install_tables():
    got = gc.genes_from_connection(install_connection(), "22")
    assert [g.stable_id for g in got] == ["G2", "G1", "G5"]  # ordered by start; other biotype and chromosome left out
    g1 = got[1]
    assert g1.exons == EXONS and g1.cds == CDS and g1.strand == 1 and g1.transcript_id == "T11" and g1.symbol == "AAA"
    assert got[0].strand == -1 and got[0].cds == () and got[0].exons == ((90, 700),)
    assert got[2].exons == () and got[2].cds == ()  # canonical transcript missing from the table
    part = gc.build_partition(got)
    assert part.utr5 == [(100, 150)] and part.utr3 == [(560, 650)]  # G2 has no CDS, so G1's UTRs stand
    assert part.exon_other == [(90, 100), (200, 300), (400, 500), (650, 700)]  # G2's non-coding exon, where nothing else claims it


def test_write_chromosome_writes_the_three_tables(tmp_path):
    sizes = gc.write_chromosome(install_connection(), "21", tmp_path / "genic_coordinates")
    assert sizes == {"cds": 30, "utr5": 5, "utr3": 10, "utr_ambiguous": 0, "exon_other": 0, "introns": 0}
    written = sorted(p.name for p in (tmp_path / "genic_coordinates").iterdir())
    assert written == ["chrom21_genes_merged.tsv", "chrom21_utr3.tsv", "chrom21_utr5.tsv"]
    header = "species\tseqid\tstart\tstop\tstrand"
    assert (tmp_path / "genic_coordinates" / "chrom21_utr5.tsv").read_text().splitlines() == [header, "homo_sapiens\t21\t5\t10\t1"]
    assert (tmp_path / "genic_coordinates" / "chrom21_utr3.tsv").read_text().splitlines() == [header, "homo_sapiens\t21\t40\t50\t1"]
    assert (tmp_path / "genic_coordinates" / "chrom21_genes_merged.tsv").read_text().splitlines() == [header, "homo_sapiens\t21\t5\t50\t1"]


def test_write_chromosome_refuses_a_chromosome_without_genes(tmp_path):
    with pytest.raises(ValueError):
        gc.write_chromosome(install_connection(), "MT", tmp_path)
    assert list(tmp_path.iterdir()) == []


def test_ambiguous_positions_are_in_neither_utr_table():
    a = gene("A", 1)  # 3' UTR (560, 650)
    c = gene("C", 1, exons=((620, 700), (800, 900)), cds=((680, 700), (800, 900)))  # 5' UTR (620, 680)
    part = gc.build_partition([a, c])
    in_tables = {p for name in ("utr5", "utr3") for _, _, start, stop, _ in gc.ref_coords_rows(part, name, "22") for p in range(start, stop)}
    assert not in_tables & positions(part.utr_ambiguous)
    assert positions(part.utr_ambiguous) == set(range(620, 650))


def test_run_writes_every_requested_chromosome(tmp_path, monkeypatch):
    monkeypatch.setattr(gc, "open_annotation", lambda install, species="human": (install_connection(), "homo_sapiens", None))
    lines = []
    written = gc.run(tmp_path, ["21", "22"], report=lines.append)
    outdir = tmp_path / "genic_coordinates"
    assert written == [outdir / f"chrom{c}_{name}.tsv" for c in ("21", "22") for name in ("utr5", "utr3", "genes_merged")]
    assert all(p.is_file() for p in written)
    assert lines[0].startswith("chromosome 21: cds 30, utr5 5, utr3 10")
    assert lines[-1].startswith("total: cds 240, utr5 55, utr3 100")
    elsewhere = gc.run(tmp_path, ["21"], outdir=tmp_path / "elsewhere")
    assert elsewhere[0] == tmp_path / "elsewhere" / "chrom21_utr5.tsv"
