from cogent3 import make_aligned_seqs, make_unaligned_seqs

from estimating_ud import apps
from tests.conftest import NAMES, cds_records, intron_records, noncds_records


def test_cds_renamer_splits_on_dash():
    seqs = make_unaligned_seqs(cds_records(), moltype="dna")
    out = apps.renamer_cds_unaligned()(seqs)
    assert sorted(out.names) == sorted(NAMES)


def test_noncds_renamer_splits_on_colon():
    aln = make_aligned_seqs(noncds_records(), moltype="dna")
    out = apps.renamer_noncds_aligned()(aln)
    assert sorted(out.names) == sorted(NAMES)


def test_wrong_renamer_drops_everything():
    # the cds renamer on ':'-named sequences finds no species prefix, so take_seqs fails
    aln = make_aligned_seqs(noncds_records(), moltype="dna")
    result = apps.renamer_cds_aligned()(aln)
    assert not result  # NotCompleted


def test_utr_samplers_cut_at_mask_boundaries():
    aln = apps.renamer_noncds_aligned()(make_aligned_seqs(intron_records(), moltype="dna"))
    utr5 = apps.sample_UTR5()(aln)
    utr3 = apps.sample_UTR3()(aln)
    body = apps.removeUTRs_fromintrons()(aln)
    assert len(utr5) == 30  # up to the first '?'
    assert len(utr3) == 300 - 260  # after the last '?'
    # aln[first '?' : last '?']: the leading masked run is kept and the trailing one dropped;
    # the '?' columns are removed later by omit_degenerates, so only the boundaries matter
    assert len(body) == 259 - 30
    assert "?" not in str(utr5.get_gapped_seq("Human"))
    assert "?" not in str(utr3.get_gapped_seq("Human"))


def test_human_seq_truncated_to_whole_motifs():
    aln = make_aligned_seqs(noncds_records(length=301), moltype="dna")
    seq = apps.gethumanseq_noncds_aligned(motif_length=3)(aln)
    assert len(seq) == 300
    # the app returns the bare sequence; stages.human_seq re-wraps it with name="Human"


def test_number_of_motifs_counts_overlapping_trinucleotides():
    from cogent3 import make_seq

    counts = apps.number_of_motifs(motif_length=3)(make_seq("ACGACG", moltype="dna", name="Human"))
    assert counts == {"ACG": 2, "CGA": 1, "GAC": 1}


def test_filter_windows_spscore_drops_gapped_window():
    good, bad = "ACGTACGTAC", "----------"
    toy = make_aligned_seqs({"Human": good * 3, "Chimpanzee": good * 3, "Orangutan": good + bad + good}, moltype="dna")
    out = apps.filter_windows_spscore(window=10, threshold_filter=0.5)(toy)
    assert out.to_dict() == {n: good * 2 for n in toy.names}
    all_gapped = make_aligned_seqs({"Human": good, "Chimpanzee": bad, "Orangutan": bad}, moltype="dna")
    assert not apps.filter_windows_spscore(window=10, threshold_filter=0.5)(all_gapped)
