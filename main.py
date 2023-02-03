import re
import sys
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
            for i, codon in enumerate(slice):
                codon_space = AMINO_ACID_LOOKUP[CODON_LOOKUP[codon]]
                out_seq.append(max(codon_space, key = lambda x: CODON_FREQUENCIES[x]))
                if codon == "CAU":
                    print(codon_space)

        return ["AUG"] + out_seq
    
    def cai(self) -> float:
        weights = [
            CODON_FREQUENCIES[x] / CODON_FREQUENCIES[max(AMINO_ACID_LOOKUP[CODON_LOOKUP[x]], key = lambda x: CODON_FREQUENCIES[x])]
            for x in self.codon_seq
        ]
        return functools.reduce(lambda x, y: x * y, weights) ** (1 / len(weights))

    @classmethod
    def from_aminoacids(cls, s: str) -> Self:
        codon_seq = [AMINO_ACID_LOOKUP[c.upper()][0] for c in s]
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


def main() -> int:
    seq = "MDIEAYLERIGYKKSRNKLDLETLTDILQHQIRAVPFENLNIHCGDAMDLGLEAIFDQVVRRNRGGWCLQVNHLLYWALTTIGFETTMLGGYVYSTPAKKYSTGMIHLLLQVTIDGRNYIVDAGFGRSYQMWQPLELISKDQPQVPCVFRLTEENGFWYLDQIRREQYIPNEEFLHSDLLEDSKYRKIYSFTLKPRTIEDFESMNTYLQTSPSSVFTSKSFCSLQTPDGVHCLVGFTLTHRRFNYKDNTDLIEFKTLSEEEIEKVLKNIFNISLQRKLVPKHGDRFFTI"
    aa_seq = AASequence.from_aminoacids(seq)
    print("".join(aa_seq.codon_seq))
    # print(aa_seq.cai())
    # print(" ".join(aa_seq.optimise()))
    
    return 0


if __name__ == "__main__":
    sys.exit(main())