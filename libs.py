import cogent3
from cogent3.app.result import model_result
from cogent3.app.composable import define_app
from phylim.apps import phylim, PhyloLimitRec

@cogent3.app.composable.define_app
def renamer(seqs: cogent3.app.typing.AlignedSeqsType) -> cogent3.app.typing.AlignedSeqsType:
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