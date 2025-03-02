#!/usr/bin/evn python3
import os
import re
from Bio import Align
from dataclasses import dataclass
from typing import Dict, Generator, Optional, Tuple

CHARGE_MAP = {
    "positive": "01",
    "negative": "10",
    None: "00",
}

@dataclass
class AminoAcid:
    name: str
    codons: set[str]
    code: str
    is_hydrophobic: bool
    is_aromatic: bool
    is_special: bool
    charge: Optional[str] = None

    def __str__(self):
        return f"{self.name} ({self.code})"

    def __repr__(self):
        return f"AminoAcid(name={self.name}, codons={self.codons}, code={self.code}, is_hydrophobic={self.is_hydrophobic}, is_aromatic={self.is_aromatic}, is_special={self.is_special}, charge={self.charge})"

    def encode_delta(self, other: "AminoAcid") -> str:
        """Encode the difference between two amino acids, treating self as mutation"""
        # Store charge
        delta_str = CHARGE_MAP[self.charge]
        delta_bits = [
            # has charge changed
            "1" if self.charge != other.charge else "0",
            # is_hydrophobic
            "1" if self.is_hydrophobic else "0",
            # has hydrophobicity changed
            "1" if other.is_hydrophobic != self.is_hydrophobic else "0",
            # is aromatic
            "1" if self.is_aromatic else "0",
            # has aromaticity changed
            "1" if other.is_aromatic != self.is_aromatic else "0",
            # is special
            "1" if self.is_special else "0",
            # has special changed
            "1" if other.is_special != self.is_special else "0",
        ]
        delta_str += "".join(delta_bits)
        return delta_str

@dataclass
class Metadata:
    title: str
    tags: Dict[str, str]

@dataclass
class Codon:
    nucleotides: str
    amino_acid: AminoAcid

AMINO_ACIDS = [
    AminoAcid(
        "Alanine",
        {"GCT", "GCC", "GCA", "GCG"},
        "A",
        is_hydrophobic=True,
        is_aromatic=False,
        is_special=False,
    ),
    AminoAcid(
        "Arginine",
        {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"},
        "R",
        charge="positive",
        is_hydrophobic=False,
        is_aromatic=False,
        is_special=False,
    ),
    AminoAcid(
        "Asparagine",
        {"AAT", "AAC"},
        "N",
        is_hydrophobic=False,
        is_aromatic=False,
        is_special=False,
    ),
    AminoAcid(
        "Aspartic Acid",
        {"GAT", "GAC"},
        "D",
        charge="negative",
        is_hydrophobic=False,
        is_aromatic=False,
        is_special=False,
    ),
    AminoAcid(
        "Cysteine",
        {"TGT", "TGC"},
        "C",
        is_hydrophobic=False,
        is_aromatic=False,
        is_special=True,
    ),
    AminoAcid(
        "Glutamine",
        {"CAA", "CAG"},
        "Q",
        is_hydrophobic=False,
        is_aromatic=False,
        is_special=False,
    ),
    AminoAcid(
        "Glutamic Acid",
        {"GAA", "GAG"},
        "E",
        charge="negative",
        is_hydrophobic=False,
        is_aromatic=False,
        is_special=False,
    ),
    AminoAcid(
        "Glycine",
        {"GGT", "GGC", "GGA", "GGG"},
        "G",
        is_hydrophobic=True,
        is_aromatic=False,
        is_special=True,
    ),
    AminoAcid(
        "Histidine",
        {"CAT", "CAC"},
        "H",
        charge="positive",
        is_hydrophobic=False,
        is_aromatic=False,
        is_special=False,
    ),
    AminoAcid(
        "Isoleucine",
        {"ATT", "ATC", "ATA"},
        "I",
        is_hydrophobic=True,
        is_aromatic=False,
        is_special=False,
    ),
    AminoAcid(
        "Leucine",
        {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"},
        "L",
        is_hydrophobic=True,
        is_aromatic=False,
        is_special=False,
    ),
    AminoAcid(
        "Lysine",
        {"AAA", "AAG"},
        "K",
        charge="positive",
        is_hydrophobic=False,
        is_aromatic=False,
        is_special=False,
    ),
    AminoAcid(
        "Methionine",
        {"ATG"},
        "M",
        is_hydrophobic=True,
        is_aromatic=False,
        is_special=False,
    ),
    AminoAcid(
        "Phenylalanine",
        {"TTT", "TTC"},
        "F",
        is_hydrophobic=True,
        is_aromatic=True,
        is_special=False,
    ),
    AminoAcid(
        "Proline",
        {"CCT", "CCC", "CCA", "CCG"},
        "P",
        is_hydrophobic=True,
        is_aromatic=False,
        is_special=False,
    ),
    AminoAcid(
        "Serine",
        {"AGT", "AGC", "TCG", "TCA", "TCC", "TCT"},
        "S",
        is_hydrophobic=False,
        is_aromatic=False,
        is_special=False,
    ),
    AminoAcid(
        "Threonine",
        {"ACT", "ACC", "ACA", "ACG"},
        "T",
        is_hydrophobic=False,
        is_aromatic=False,
        is_special=False,
    ),
    AminoAcid(
        "Tryptophan",
        {"TGG"},
        "W",
        is_hydrophobic=True,
        is_aromatic=True,
        is_special=False,
    ),
    AminoAcid(
        "Tyrosine",
        {"TAT", "TAC"},
        "Y",
        is_hydrophobic=True,
        is_aromatic=True,
        is_special=False,
    ),
    AminoAcid(
        "Valine",
        {"GTT", "GTC", "GTA", "GTG"},
        "V",
        is_hydrophobic=True,
        is_aromatic=False,
        is_special=False,
    ),
]

CODON_TO_AMINO_ACID = {
    codon: amino_acid
    for amino_acid in AMINO_ACIDS
    for codon in amino_acid.codons
}

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(THIS_DIR, "test_data")

TEST_DATA_PATHS = {
    "ecoli_isolated_protein_ref": os.path.join(DATA_DIR, "ncbi_lacZ_dataset/data/gene.fna"),
    "ecoli_isolated_protein_alpha": os.path.join(DATA_DIR, "ncbi_lacZ_dataset/data/gene-5.fna"),
    "ecoli_isolated_protein_beta": os.path.join(DATA_DIR, "ncbi_lacZ_dataset/data/gene-6.fna"),
}

def load_fasfa_samples(path: str) -> Generator[Tuple[Metadata, str], None, None]:
    """Load fasta samples from a file"""
    with open(path, "r") as f:
        metadata = None
        sequence = ""
        for line in f:
            if line.startswith(">"):
                if metadata is not None:
                    yield metadata, sequence

                key_value_pairs = re.findall(r"\[(.*?)=(.*?)\]", line)
                tags = {k: v for k, v in key_value_pairs}

                metadata = Metadata(
                    title=line[1:line.index("[")].strip(),
                    tags=tags,
                )

                sequence = ""
            else:
                sequence += line.strip()

        if sequence:
            yield metadata, sequence

def chunk_protein(sequence: str) -> Generator[str, None, None]:
    """Given a protein starting in ATG and ending in TAA, TGA, or TAG, yield chunks of 3"""
    if not sequence.startswith("ATG"):
        raise ValueError("Sequence must start with ATG")

    if not sequence.endswith(("TAA", "TGA", "TAG")):
        raise ValueError("Sequence must end with TAA, TGA, or TAG")

    for i in range(0, len(sequence), 3):
        yield sequence[i:i+3]

def is_start_codon(codon: str) -> bool:
    """Given a codon, return True if it is a start codon"""
    return codon == "ATG"

def is_stop_codon(codon: str) -> bool:
    """Given a codon, return True if it is a stop codon"""
    return codon in ("TAA", "TGA", "TAG")

def map_codon(codon: str) -> Optional[AminoAcid]:
    """Maps a codon string to an amino acid"""
    return CODON_TO_AMINO_ACID.get(codon)

def main():
    ecoli_ref_data = TEST_DATA_PATHS["ecoli_isolated_protein_ref"]
    ecoli_sample_data = TEST_DATA_PATHS["ecoli_isolated_protein_alpha"]

    ref_sequence = next(load_fasfa_samples(ecoli_ref_data))
    sample_sequence = next(load_fasfa_samples(ecoli_sample_data))

    print(f"Reference sequence: {ref_sequence[0].title}, {ref_sequence[0].tags}")
    print(f"Sample sequence: {sample_sequence[0].title}, {sample_sequence[0].tags}")

    aligner = Align.PairwiseAligner(scoring="blastn")
    alignments = aligner.align(ref_sequence[1], sample_sequence[1])

    # Print to compare sequence lengths
    print(f"Reference sequence length: {len(ref_sequence[1])}")
    print(f"Sample sequence length: {len(sample_sequence[1])}")

    for alignment in alignments:
        print(f"Score: {alignment.score}")

        # Print aligned segments of 50 chars from each on top of the other
        # until we have printed the whole of each sequence

        print(f"Aligned segments:")
        print(f"{alignment[0]}")
        print(f"{alignment[1]}")

        break

    print("Amino acid encoding:")
    print(f"Amino acid {AMINO_ACIDS[11]} mutating to {AMINO_ACIDS[9]}")
    print(AMINO_ACIDS[9].encode_delta(AMINO_ACIDS[11]))

if __name__ == "__main__":
    main()
