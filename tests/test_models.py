import importlib

from cogent3 import get_app

from estimating_ud import dinuc_models, models

GN = ["A>C", "A>T", "A>G", "C>A", "C>T", "C>G", "T>A", "T>C", "G>A", "G>C", "G>T"]


def test_trinuc_models_registered_and_param_lists_unchanged():
    # guards against the GN/CpG helper functions being swapped for the dinucleotide ones
    assert models.GT_CpG_ss().get_param_list() == GN + ["(CG>TG | CG>CA)", "omega"]
    assert models.GNC_CpG_ss().get_param_list() == GN + ["(CG>TG | CG>CA)", "omega"]
    assert get_app("model", "GT_CpG_ss", show_progress=False)
    assert get_app("model", "GNC_CpG_ss", show_progress=False)


def test_dinuc_models_param_lists():
    assert dinuc_models.GDN_CpG_ss().get_param_list() == GN + ["(CG>TG | CG>CA)"]
    assert dinuc_models.GDN_CpG().get_param_list() == GN + ["CG>TG", "CG>CA"]


def test_reimport_does_not_raise():
    importlib.reload(models)
    assert get_app("model", "GT_CpG_ss", show_progress=False)


def test_stage5_alphabet_matches_model_motifs():
    from estimating_ud.stages.constraint import trinuc_alphabet

    alphabet = list(trinuc_alphabet())
    assert len(alphabet) == 64
    assert set(models.GT_CpG_ss().get_motifs()) == set(alphabet)
    # the codon model lacks exactly the three stop codons that stage 5 pads in
    assert set(alphabet) - set(models.GNC_CpG_ss().get_motifs()) == {"TAA", "TAG", "TGA"}
