"""End-to-end run of the stages on the synthetic data root (chromosome 22 only)."""

import pickle

import pytest
from cogent3 import load_aligned_seqs, load_unaligned_seqs

from estimating_ud import regions
from estimating_ud.stages import align_cds, constraint, filter_noncds, fit, human_seq, quality


def alnstat(path):
    return constraint.parse_alignment_length(path.read_text())


@pytest.mark.parametrize("submodel", ["singlent", "trinuc"])
def test_align_cds_writes_filtered_alignment(data_root, submodel):
    out = align_cds.run("22", submodel, data_root)
    aln = load_aligned_seqs(str(out), moltype="dna")
    assert sorted(aln.names) == ["Chimpanzee", "Human", "Orangutan"]
    assert len(aln) % 3 == 0
    assert len(aln) > 0
    assert alnstat(regions.alnstat_path(data_root, "cds", "22", submodel)) == f"{len(aln):.6e}"


def test_filter_noncds_drops_masked_columns_of_introns(data_root):
    # introns arrive with exons and repeats already masked; stage 3 only has to drop those columns
    (out,) = filter_noncds.run("22", "singlent", data_root, regions=("introns",))
    aln = load_aligned_seqs(str(out), moltype="dna")
    text = "".join(str(s) for s in aln.seqs)
    assert "?" not in text and "-" not in text
    assert len(aln) == 2 * (300 - 20 - 6)  # two fragments, each with 20 masked and 6 gapped columns


def test_filter_noncds_trinuc_introns_keeps_whole_motifs(data_root):
    (out,) = filter_noncds.run("22", "trinuc", data_root, regions=("introns",))
    aln = load_aligned_seqs(str(out), moltype="dna")
    assert len(aln) > 0 and len(aln) % 3 == 0


def test_gene_named_intron_fragments_are_refused(tmp_path):
    # the per-gene download used before introns were defined by coordinates must not be filtered as the new class
    from tests.conftest import intron_records, write_fasta

    write_fasta(tmp_path / "introns" / "alldata_chrm22" / "ENSG00000100001.fa", intron_records())
    with pytest.raises(ValueError, match="not named by coordinates"):
        filter_noncds.run("22", "singlent", tmp_path, regions=("introns",))
    with pytest.raises(ValueError, match="not named by coordinates"):
        human_seq.run("22", tmp_path, regions=("introns",))


def test_filter_noncds_two_regions(data_root):
    outs = filter_noncds.run("22", "singlent", data_root, regions=("intergenicAR", "distalIG"))
    for p in outs:
        aln = load_aligned_seqs(str(p), moltype="dna")
        assert "-" not in "".join(str(s) for s in aln.seqs)
        assert 0 < len(aln) <= 600


def test_human_seq_concatenates(data_root):
    outs = human_seq.run("22", data_root, regions=("cds", "introns", "intergenicAR"))
    for p in outs:
        seq = load_unaligned_seqs(str(p), moltype="dna").seqs["Human"]
        assert len(seq) > 0
        assert "?" not in str(seq) or p.parent.parent.name == "cds"


def test_fit_and_constraint_singlent(data_root):
    filter_noncds.run("22", "singlent", data_root, regions=("intergenicAR", "distalIG"))
    align_cds.run("22", "singlent", data_root)
    pickles = fit.run("22", "singlent", data_root, regions=("intergenicAR", "distalIG", "cds"))
    for p in pickles:
        with open(p, "rb") as f:
            result = pickle.load(f)
        assert result.lf.get_rate_matrix_for_edge("Human", calibrated=False).shape == (4, 4)

    ens, cons = constraint.compute(data_root, "singlent", chromosomes=("22",), regions=("distalIG", "cds"))
    assert list(ens["Region"]) == ["IGAR", "distalIG", "distalIG", "cds", "cds"]
    assert list(cons["Region"]) == ["distalIG", "cds"]
    assert cons["Constraint"].notna().all()


def test_fit_trinuc_cds_pads_stop_codons(data_root):
    align_cds.run("22", "trinuc", data_root)
    filter_noncds.run("22", "trinuc", data_root, regions=("intergenicAR",))
    fit.run("22", "trinuc", data_root, regions=("intergenicAR", "cds"))
    ens, cons = constraint.compute(data_root, "trinuc", chromosomes=("22",), regions=("cds",))
    assert ens.shape == (3, 5)
    assert cons.shape == (1, 3)


def test_quality_branch_ens(data_root, tmp_path):
    # the fixture has no intronsAR, so only the region the later steps read
    quality.masked_databank("22", data_root, regions=("intergenicAR",))
    (fasta,) = quality.filter_spscore("22", data_root)
    quality.fit_spscore("22", data_root)
    aln = load_aligned_seqs(str(fasta), moltype="dna")
    assert len(aln) > 0 and len(aln) % 3 == 0

    ens = quality.ens_spscore(data_root, chromosomes=("22",))
    assert list(ens.columns) == quality.ENS_COLUMNS
    assert ens.shape == (1, 6)
    row = ens.iloc[0]
    assert (row["Region"], row["Chromosome"]) == ("IGAR", "22")
    assert (row["window"], row["threshold"]) == (regions.QUALITY_WINDOW, regions.QUALITY_THRESHOLD)
    assert row["aln_length"] == len(aln)
    assert 0 < row["ENS"] < 1

    out = quality.write_ens_spscore(data_root, tmp_path / "tables", chromosomes=("22",))
    assert out.name == "trinuc_ENS_spscore_window60_qualitythreshold0.5.csv"
    assert out.read_text().splitlines()[0] == ",".join(quality.ENS_COLUMNS)

    with pytest.raises(FileNotFoundError):
        quality.ens_spscore(data_root, window=30, chromosomes=("22",))
