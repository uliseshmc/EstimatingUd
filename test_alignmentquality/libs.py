import cogent3
from cogent3.app.result import model_result
from cogent3.app.composable import define_app
from phylim.apps import phylim, PhyloLimitRec

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

@cogent3.app.composable.define_app
def degap_app(seqs: cogent3.app.typing.AlignedSeqsType) -> cogent3.app.typing.SeqsCollectionType:
    """
    A function to degap aligned sequences and return an unaligned sequence.
    """

    return seqs.degap()

