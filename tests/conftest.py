"""Synthetic fixtures: a fake data root with one tiny fragment per region family.

Sequence naming mirrors eti output: cds homologs are ``<species>-<gene id>`` (unaligned),
every other region is ``<species>:<chrom>:<start>-<end>:<strand>`` (aligned). Intron
fragments carry ``?`` runs in the Human sequence where exons / repeats were masked.
"""

from __future__ import annotations

import random
from pathlib import Path

import pytest

SPECIES = ("homo_sapiens", "pan_troglodytes", "pongo_abelii")
NAMES = ("Human", "Chimpanzee", "Orangutan")


def _mutate(seq: str, n: int, seed: int) -> str:
    rng = random.Random(seed)
    s = list(seq)
    for _ in range(n):
        i = rng.randrange(len(s))
        s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


def cds_records(n_codons: int = 40, seed: int = 1) -> dict[str, str]:
    """Three unaligned, stop-terminated coding sequences without internal stops."""
    rng = random.Random(seed)
    sense = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in ("TAA", "TAG", "TGA")]
    human = "".join(rng.choice(sense) for _ in range(n_codons))
    seqs = {}
    for i, sp in enumerate(SPECIES):
        s = human
        # a few synonymous-position changes, avoiding stop creation by re-checking codons
        while True:
            cand = _mutate(s, 2 + 3 * i, seed + i)
            if all(cand[j : j + 3] not in ("TAA", "TAG", "TGA") for j in range(0, len(cand), 3)):
                break
            seed += 10
        seqs[f"{sp}-ENSG{i:011d}"] = cand + "TAA"
    return seqs


def noncds_records(length: int = 300, seed: int = 2, chrom: str = "22", start: int = 1000) -> dict[str, str]:
    rng = random.Random(seed)
    human = "".join(rng.choice("ACGT") for _ in range(length))
    seqs = {}
    for i, sp in enumerate(SPECIES):
        s = _mutate(human, 3 + 4 * i, seed + i)
        if i == 2:  # one gap run in the outgroup so omit_degenerates has work to do
            s = s[:60] + "-" * 6 + s[66:]
        seqs[f"{sp}:{chrom}:{start}-{start + length}:1"] = s
    return seqs


def intron_records(length: int = 300, seed: int = 3) -> dict[str, str]:
    """Like noncds_records, with masked (?) runs at columns 30-39 and 250-259 in the Human sequence."""
    seqs = noncds_records(length, seed)
    human_key = next(k for k in seqs if k.startswith("homo_sapiens"))
    h = list(seqs[human_key])
    for a, b in ((30, 40), (250, 260)):
        h[a:b] = "?" * (b - a)
    seqs[human_key] = "".join(h)
    return seqs


def write_fasta(path: Path, records: dict[str, str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records.items()))
    return path


@pytest.fixture(scope="session")
def data_root(tmp_path_factory) -> Path:
    """A fake alignment root with raw fragments for cds, introns and intergenicAR on chr22."""
    root = tmp_path_factory.mktemp("data_root")
    chrm = "22"
    for i in range(2):
        write_fasta(root / "cds" / f"alldata_chrm{chrm}" / f"gene{i}.fa", cds_records(seed=10 + i))
        # coordinate-defined regions are named by position, as eti alignments --ref_coords writes them
        write_fasta(
            root / "introns" / f"alldata_chrm{chrm}" / f"homo_sapiens-{chrm}-{1000 * i}-{1000 * i + 300}.fa",
            intron_records(seed=20 + i),
        )
        write_fasta(root / "intergenicAR" / f"alldata_chrm{chrm}" / f"block{i}.fa", noncds_records(seed=30 + i, start=1000 * i))
        write_fasta(root / "distalIG" / f"alldata_chrm{chrm}" / f"block{i}.fa", noncds_records(seed=40 + i, start=1000 * i))
    return root
