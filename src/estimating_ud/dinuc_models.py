"""Dinucleotide substitution models with CpG deamination (``GDN_CpG``, ``GDN_CpG_ss``).

Kept apart from :mod:`estimating_ud.models` on purpose: both modules define helpers named
``make_gn_preds`` etc., and their bodies differ (this one enumerates 11 ``MotifChange``
predicates, the trinucleotide module uses cogent3's ``_general_preds``). Merging them would
silently change the trinucleotide models.
"""

from itertools import permutations

from cogent3.evolve.ns_substitution_model import NonReversibleDinucleotide
from cogent3.evolve.predicate import MotifChange


def make_gn_preds():
    # making the model parameters (predictates) for
    # the General Nucleotide Markov model
    return [
        MotifChange(f, t, forward_only=True)
        for f, t in permutations("ACTG", 2)
        if f != "T" or t != "G"
    ]


def make_nr_cpg_preds_strand_asymetric():
    # strand specific CpG deamination rates, so separate
    # parameters for plus strand (CG->TG) and minus
    # strand (CG->CA)
    return [
        MotifChange("CG", "TG", forward_only=True),
        MotifChange("CG", "CA", forward_only=True),
    ]


def make_nr_cpg_preds_strand_symetric():
    # same CpG deamination rate on both strands
    # so one parameter
    return [
        # | is the binary or operator that combines the predicates
        # so we have the union of the two changes as a single parameter
        MotifChange("CG", "TG", forward_only=True) | MotifChange("CG", "CA", forward_only=True)
    ]


def _make_model(cls, **kwargs):
    return cls(**kwargs)


def GDN_CpG_ss(**kwargs):
    """return a dinucleotide model with strand symmetric CpG deamination"""
    ssym_preds = make_gn_preds() + make_nr_cpg_preds_strand_symetric()
    kwargs = dict(predicates=ssym_preds, optimise_motif_probs=True, name="GDN_CpG_ss")
    return _make_model(NonReversibleDinucleotide, **kwargs)


def GDN_CpG(**kwargs):
    """return a dinucleotide model with strand asymmetric CpG deamination"""
    asym_preds = make_gn_preds() + make_nr_cpg_preds_strand_asymetric()
    kwargs = dict(predicates=asym_preds, optimise_motif_probs=True, name="GDN_CpG")
    return _make_model(NonReversibleDinucleotide, **kwargs)
