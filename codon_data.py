# import itertools

# CODONS = ["".join(x) for x in itertools.product("AUGC", "AUGC", "AUGC")]

START_CODONS = {"AUG"}

STOP_CODONS = {"UAA", "UAG", "UGA"}

AMINO_ACID_NAMES = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Asx": "B",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Glx": "Z",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V"
}

AMINO_ACID_LOOKUP = {
    "A": ["GCU", "GCC", "GCA", "GCG"],
    "C": ["UGU", "UGC"],
    "D": ["GAU", "GAC"],
    "E": ["GAA", "GAG"],
    "F": ["UUU", "UUC"],
    "G": ["GGU", "GGC", "GGA", "GGG"],
    "H": ["CAU", "CAC"],
    "I": ["AUU", "AUC", "AUA"],
    "K": ["AAA", "AAG"],
    "L": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
    "M": ["AUG"],
    "N": ["AAU", "AAC"],
    "P": ["CCU", "CCC", "CCA", "CCG"],
    "Q": ["CAA", "CAG"],
    "R": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "S": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
    "T": ["ACU", "ACC", "ACA", "ACG"],
    "V": ["GUU", "GUC", "GUA", "GUG"],
    "W": ["UGG"],
    "Y": ["UAU", "UAC"]
}

CODON_LOOKUP = {}
for k, v in AMINO_ACID_LOOKUP.items():
    CODON_LOOKUP.update(
        {codon: k for codon in v}
    )

CODON_FREQUENCIES = {
    "AAA": 0.43,
	"AAC": 0.53,
	"AAG": 0.57,
	"AAU": 0.47,
	"ACA": 0.28,
	"ACC": 0.36,
	"ACG": 0.11,
	"ACU": 0.25,
	"AGA": 0.21,
	"AGC": 0.24,
	"AGG": 0.21,
	"AGU": 0.15,
	"AUA": 0.17,
	"AUC": 0.47,
	"AUG": 1.0,
	"AUU": 0.36,
	"CAA": 0.27,
	"CAC": 0.58,
	"CAG": 0.73,
	"CAU": 0.42,
	"CCA": 0.28,
	"CCC": 0.32,
	"CCG": 0.11,
	"CCU": 0.29,
	"CGA": 0.11,
	"CGC": 0.18,
	"CGG": 0.2,
	"CGU": 0.08,
	"CUA": 0.07,
	"CUC": 0.2,
	"CUG": 0.4,
	"CUU": 0.13,
	"GAA": 0.42,
	"GAC": 0.54,
	"GAG": 0.58,
	"GAU": 0.46,
	"GCA": 0.23,
	"GCC": 0.4,
	"GCG": 0.11,
	"GCU": 0.27,
	"GGA": 0.25,
	"GGC": 0.34,
	"GGG": 0.25,
	"GGU": 0.16,
	"GUA": 0.12,
	"GUC": 0.24,
	"GUG": 0.46,
	"GUU": 0.18,
	"UAA": 0.3,
	"UAC": 0.56,
	"UAG": 0.24,
	"UAU": 0.44,
	"UCA": 0.15,
	"UCC": 0.22,
	"UCG": 0.05,
	"UCU": 0.19,
	"UGA": 0.47,
	"UGC": 0.54,
	"UGG": 1.0,
	"UGU": 0.46,
	"UUA": 0.08,
	"UUC": 0.54,
	"UUG": 0.13,
	"UUU": 0.46
}

