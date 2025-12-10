import cogent3
from cogent3.app.result import model_result
from cogent3.app.composable import define_app
from phylim.apps import phylim, PhyloLimitRec

from itertools import permutations
import cogent3 as c3
from cogent3.evolve.predicate import MotifChange
from cogent3.evolve.ns_substitution_model import NonReversibleDinucleotide, NonReversibleCodon

# I'm using this app to rename sequences in a dataset
# Some alignment files include duplicates (paralogs). This function will throw an error and not raneme the sequences if it happens

@cogent3.app.composable.define_app
def renamer_cds_unaligned(seqs: cogent3.app.typing.UnalignedSeqsType) -> cogent3.app.typing.UnalignedSeqsType:
    """
    A function to rename sequences in a dataset.
    """
    name_map = {
        "homo_sapiens": "Human",
        "pan_troglodytes": "Chimpanzee",
        "gorilla_gorilla": "Gorilla"
    }

    seqs = seqs.rename_seqs(lambda x: name_map.get(x.split("-")[0], x))

    return seqs.take_seqs(list(name_map.values()))

@cogent3.app.composable.define_app
def renamer_cds_aligned(seqs: cogent3.app.typing.AlignedSeqsType) -> cogent3.app.typing.AlignedSeqsType:
    """
    A function to rename sequences in a dataset.
    """
    name_map = {
        "homo_sapiens": "Human",
        "pan_troglodytes": "Chimpanzee",
        "gorilla_gorilla": "Gorilla"
    }

    seqs = seqs.rename_seqs(lambda x: name_map.get(x.split("-")[0], x))

    return seqs.take_seqs(list(name_map.values()))

@cogent3.app.composable.define_app
def renamer_noncds_unaligned(seqs: cogent3.app.typing.UnalignedSeqsType) -> cogent3.app.typing.UnalignedSeqsType:
    """
    A function to rename sequences in a dataset.
    """
    name_map = {
        "homo_sapiens": "Human",
        "pan_troglodytes": "Chimpanzee",
        "gorilla_gorilla": "Gorilla"
    }

    seqs = seqs.rename_seqs(lambda x: name_map.get(x.split(":")[0], x))

    return seqs.take_seqs(list(name_map.values()))

@cogent3.app.composable.define_app
def renamer_noncds_aligned(seqs: cogent3.app.typing.AlignedSeqsType) -> cogent3.app.typing.AlignedSeqsType:
    """
    A function to rename sequences in a dataset.
    """
    name_map = {
        "homo_sapiens": "Human",
        "pan_troglodytes": "Chimpanzee",
        "gorilla_gorilla": "Gorilla"
    }

    seqs = seqs.rename_seqs(lambda x: name_map.get(x.split(":")[0], x))

    return seqs.take_seqs(list(name_map.values()))

#this is a workaround to use phylim on models with splitted codons
#I'm using a workaround to check identifiability of a nucleotide model split by position
#Latter I should check this workaround in case I actually use this model
#Phylim is getting updated to fix such bug
@define_app
def phylim_split_codon(result: model_result, check_one: phylim) -> PhyloLimitRec:
    """checks individual likelihood functions from a split_codon model_result"""
    for k in range (1, 4):
        value = result[k]
        one = model_result(name=result.name, source=result.source)
        one['value'] = value
        checked = check_one(one)
        if not checked.is_identifiable:
            return checked
    return checked

#Next functions define a dinucleotide substitution model. This model is useful for CpG evolution sites. 
#The way to call the model is by using 
# subsmodel = c3.get_app("model", GDN_CpG_ss(), time_het="max", show_progress=True, optimise_motif_probs=True)
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
        MotifChange("CG", "TG", forward_only=True) |
        MotifChange("CG", "CA", forward_only=True)
    ]

def _make_model(cls, **kwargs):
    return cls(**kwargs)

def GDN_CpG_ss(**kwargs):
    """return a dinucleotide model with strand symmetric CpG deamination"""
    ssym_preds = make_gn_preds() + make_nr_cpg_preds_strand_symetric()
    kwargs=dict(predicates=ssym_preds, optimise_motif_probs=True, name="GDN_CpG_ss")
    return _make_model(NonReversibleDinucleotide, **kwargs)

def GDN_CpG(**kwargs):
    """return a dinucleotide model with strand asymmetric CpG deamination"""
    asym_preds = make_gn_preds() + make_nr_cpg_preds_strand_asymetric()
    kwargs=dict(predicates=asym_preds, optimise_motif_probs=True, name="GDN_CpG")
    return _make_model(NonReversibleDinucleotide, **kwargs)