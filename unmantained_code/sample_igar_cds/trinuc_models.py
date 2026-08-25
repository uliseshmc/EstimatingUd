from cogent3.evolve.models import _general_preds, register_model
from cogent3.evolve.predicate import MotifChange, UserPredicate
from cogent3.evolve.ns_substitution_model import (
    NonReversibleCodon,
    NonReversibleTrinucleotide,
)

# to build omega for both trinuc and codon
from cogent3 import get_code
from cogent3.evolve.substitution_model import _CodonPredicates
import typing


if typing.TYPE_CHECKING:
    from cogent3.core.genetic_code import GeneticCode


def make_omega_preds():
    # making the model parameters (predictates) for
    # the omega parameter of codon models
    gc = get_code(1)
    codon_preds = _CodonPredicates(gc)
    return UserPredicate(codon_preds.replacement).aliased("omega")


omega = [make_omega_preds()]


def make_gn_preds():
    # making the model parameters (predictates) for
    # the General Nucleotide Markov model
    return _general_preds


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
        MotifChange("CG", "TG", forward_only=True)
        | MotifChange("CG", "CA", forward_only=True)
    ]


def _make_model(cls, **kwargs):
    return cls(**kwargs)


@register_model("codon")
def GNC_CpG_ss(**kwargs):
    """return a CODON model with GN predicates, omega, and strand symmetric CpG deamination"""
    ssym_preds = make_gn_preds() + make_nr_cpg_preds_strand_symetric() + omega
    cpgdecay = ssym_preds[-1]
    kwargs = dict(
        predicates=ssym_preds,
        optimise_motif_probs=True,
        name="GNC_CpG_ss",
        scales={"CpG": cpgdecay, "notCpG": ~cpgdecay},
    )
    return _make_model(NonReversibleCodon, **kwargs)


@register_model("codon")
def GT_CpG_ss(**kwargs):
    """return a Trinucleotide model with GN predicates, omega, and strand
    symmetric CpG deamination

    Notes
    -----
    This is identical to GNC_CpG_ss except that it has
    ALL possible trinucleotides as states, rather than just the
    sense codons.
    """
    ssym_preds = make_gn_preds() + make_nr_cpg_preds_strand_symetric() + omega
    #print(ssym_preds)
    cpgdecay = ssym_preds[-1]
    kwargs = dict(
        predicates=ssym_preds,
        optimise_motif_probs=True,
        name="GT_CpG_ss",
        scales={"CpG": cpgdecay, "notCpG": ~cpgdecay},
    )
    return _make_model(NonReversibleTrinucleotide, **kwargs)


def modified_lf(lf, gc: int | str | None = None):
    """return a new likelihood function with the fitted parameter rules

    Parameters
    ----------
    lf
        a likelihood function that has been optimised
    """
    rules = lf.get_param_rules()
    new_lf = lf.model.make_likelihood_function(tree=lf.tree)
    new_lf.apply_param_rules(rules)
    return new_lf
