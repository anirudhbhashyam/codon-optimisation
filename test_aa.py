from typing import NoReturn

import pytest

from main import *


CASES = [
    (
        "MVLHEQSGSTLL",
        "AUGGUGCUGCACGAGCAGAGCGGCAGCACCCUGCUG"
    ),
    (
        "MDIEAYLERIGYKKSRNKLDLETLTDILQHQIRAVPFENLNIHCGDAMDLGLEAIFDQVVRRNRGGWCLQVNHLLYWALTTIGFETTMLGGYVYSTPAKKYSTGMIHLLLQVTIDGRNYIVDAGFGRSYQMWQPLELISGKDQPQVPCVFRLTEENGFWYLDQIRREQYIPNEEFLHSDLLEDSKYRKIYSFTLKPRTIEDFESMNTYLQTSPSSVFTSKSFCSLQTPDGVHCLVGFTLTHRRFNYKDNTDLIEFKTLSEEEIEKVLKNIFNISLQRKLVPKHGDRFFTI",
        ""
    )
]


@pytest.mark.parametrize(
    ("seq", "expected"), 
    [
        ("ACDEFGHIKLMNPQRSTVWY", "GCUUGUGAUGAAUUUGGUCAUAUUAAAUUAAUGAAUCCUCAACGUUCUACUGUUUGGUAU"),
        (
            "MVLHEQSGSTLLVENTESAIEINDSQSSKDPVRLALLYSEDMGAGCDINRHAESTLVGLSNKTANANCAMIMHPVVVDLPVRNADSRMSEMCLSKIKGEP", 
            "AUGGUUUUACAUGAACAAUCUGGUUCUACUUUAUUAGUUGAAAAUACUGAAUCUGCUAUUGAAAUUAAUGAUUCUCAAUCUUCUAAAGAUCCUGUUCGUUUAGCUUUAUUAUAUUCUGAAGAUAUGGGUGCUGGUUGUGAUAUUAAUCGUCAUGCUGAAUCUACUUUAGUUGGUUUAUCUAAUAAAACUGCUAAUGCUAAUUGUGCUAUGAUUAUGCAUCCUGUUGUUGUUGAUUUACCUGUUCGUAAUGCUGAUUCUCGUAUGUCUGAAAUGUGUUUAUCUAAAAUUAAAGGUGAACCU"
        )
    ]
)
def test_conversion(seq: str, expected: str) -> NoReturn:
    aa_seq = AASequence.from_aminoacids(seq)
    assert "".join(aa_seq.codon_seq) == expected


@pytest.mark.parametrize(
    ("seq", "expected"),
    CASES
)
def test_aa(seq: str, expected: str) -> NoReturn:
    aa_seq = AASequence.from_aminoacids(seq)
    assert "".join(aa_seq.optimise()) == expected
    