import re
import sys
import argparse
import itertools
import functools

from dataclasses import dataclass

from typing import Generator
from typing import Iterator

try:
    from typing import Self
except ImportError:
    from typing import TypeVar
    Self = TypeVar("Self")

from codon_data import *

take_tuples = lambda g, n = 3: zip(*(itertools.islice(gen, i, None, n) for i, gen in enumerate(itertools.tee(g, n))))
# start_stop_regex = r"(?<=(AUG)).+(?=[(UAA)|(UAG)|(UGA)]"

@dataclass
class AASequence:
    codon_seq: Iterator[str]
    start_pos: Generator[int, None, None]

    def optimise(self) -> str:
        out_seq = []
        for slice in self._seq_iter():
            new_seq = []
            for codon in slice:
                codon_space = AMINO_ACID_LOOKUP[CODON_LOOKUP[codon]]
                new_seq.append(max(codon_space, key = lambda x: CODON_FREQUENCIES[x]))
            out_seq.append(["AUG"] + new_seq)
        return out_seq
    
    def cai(self) -> float:
        weights = [
            CODON_FREQUENCIES[x] / CODON_FREQUENCIES[max(AMINO_ACID_LOOKUP[CODON_LOOKUP[x]], key = lambda x: CODON_FREQUENCIES[x])]
            for x in self.codon_seq
        ]
        return functools.reduce(lambda x, y: x * y, weights) ** (1 / len(weights))

    @classmethod
    def from_aminoacids_A(cls, s: str) -> Self:
        codon_seq = [AMINO_ACID_LOOKUP[c.upper()][0] for c in s]
        start_indices = (i for i, x in enumerate(codon_seq) if x in START_CODONS)
        return cls(codon_seq, start_indices)
    
    @classmethod
    def from_aminoacids_a(cls, s: str) -> Self:
        aa_seq = (AMINO_ACID_NAMES[name] for name in take_tuples(s))
        codon_seq = [AMINO_ACID_LOOKUP[l][0] for l in aa_seq]
        start_indices = (i for i, x in enumerate(codon_seq) if x in START_CODONS)
        return cls(codon_seq, start_indices)
    
    @classmethod
    def from_codons(cls, s: str, frame: int = None) -> Self:
        if frame is None:
            frame = 0
        start_indices = (x.start() for x in re.finditer(r"(AUG)", s[frame :]))
        codon_seq = (
            ("".join(triplet) for triplet in itertools.takewhile(lambda x: x not in STOP_CODONS, take_tuples(itertools.islice(s[frame :], start_index, None)))) 
            for start_index in start_indices
        )
        return cls(itertools.chain(*codon_seq), start_indices)
    
    def _seq_iter(self) -> Generator[Generator[str, None, None], None, None]:
        for start_index in self.start_pos:
            seq_slice = itertools.islice(self.codon_seq, start_index + 1, None)
            yield (triplet for triplet in itertools.takewhile(lambda x: x not in STOP_CODONS & START_CODONS, seq_slice))


def process_args() -> argparse.Namespace:
    args = argparse.ArgumentParser()
    args.add_argument(
        "mode",
        type = str,
        help = "The type of sequence provided."
    )
    args.add_argument(
        "sequence",
        type = str,
        help = "The sequence to optimise."
    )
    return args.parse_args()


def main() -> int:
    args = process_args()
    match args.mode:
        case "A":
            aa_seq = AASequence.from_aminoacids_A(args.sequence)
        case "a":
            aa_seq = AASequence.from_aminoacids_a(args.sequence)
        case _:
            raise ValueError(f"Sequence type {args.mode} not supported. Please give an amino acid sequence in reduced form.")
    # print("".join(aa_seq.codon_seq))
    print(["".join(slice) for slice in aa_seq.optimise()])
    return 0


if __name__ == "__main__":
    sys.exit(main())