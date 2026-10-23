"""End-to-end run of the stages on the synthetic data root (chromosome 22 only)."""

import pickle

import pytest
from cogent3 import load_aligned_seqs, load_unaligned_seqs

from estimating_ud import regions
from estimating_ud.stages import align_cds, constraint, filter_introns, filter_noncds, fit, human_seq


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


def test_filter_introns_samples_three_regions(data_root):
    outs = filter_introns.run("22", "trinuc", data_root)
    lengths = {p.parent.parent.name: len(load_aligned_seqs(str(p), moltype="dna")) for p in outs}
    assert set(lengths) == set(regions.INTRON_REGIONS)
    assert all(v > 0 and v % 3 == 0 for v in lengths.values())
    assert lengths["introns_nonUTR"] > lengths["introns5UTR"]


def test_filter_noncds_two_regions(data_root):
    outs = filter_noncds.run("22", "singlent", data_root, regions=("intergenicAR", "distalIG"))
    for p in outs:
        aln = load_aligned_seqs(str(p), moltype="dna")
        assert "-" not in "".join(str(s) for s in aln.seqs)
        assert 0 < len(aln) <= 600


def test_human_seq_concatenates(data_root):
    outs = human_seq.run("22", data_root, regions=("cds", "introns_nonUTR", "intergenicAR"))
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
