# Codon Optimisation

A simple algorithm to implement codon optimisation for humans.

## CI

![test](https://www.github.com/anirudhbhashyam/codon-optimisation/actions/workflows/test.yml/badge.svg)

# Usage
```sh
$ git clone "https://www.github.com/anirudhbhashyam/codon-optimisation"
$ cd codon-optimisation
$ python main.py <sequence-type> <sequence>
```

## Arugments

- `sequence-type`: Supported types are `A` and `a`. `A` represents a sequence of amino acids provided using name abbreviations (eg: "MVB..."). `a` represents a sequence of amino acids supplied in short form (eg: "TyrMetCys...")
- `sequence`: A sequence of amino acids.

# Dependencies

## Use
- `python >= 3.10`

## Development
- `pytest`

# References

- Codon Usage Frequency: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606&aa=1&style=N
- Random Amino Acid Sequences: https://molbiotools.com/randomsequencegenerator.php


