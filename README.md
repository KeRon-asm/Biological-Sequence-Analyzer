# FASTA Sequence Analyzer

A bioinformatics tool written in C++ for parsing, analyzing, and translating biological sequences from FASTA files.

## Features

- **FASTA Parsing** — Reads and extracts sequences and headers from `.fasta` files
- **Sequence Statistics** — Computes length, amino acid composition, and molecular weight
- **Motif Search** — Searches for sequence motifs using a Trie data structure
- **DNA Translation** — Translates DNA sequences to protein sequences using a full codon table
- **Local Alignment** — Smith-Waterman local alignment algorithm *(in progress)*

## Project Structure

```
fasta-analyzer/
├── main.cpp        # Entry point, coordinates program flow
├── fasta.h         # Declarations for structs and functions
├── fasta.cpp       # Implementations for parsing and codon table
└── EXAMPLE_1.fasta # Sample FASTA file for testing
```

## Getting Started

### Prerequisites

- g++ with C++17 support or later
- A terminal (Linux, macOS, or WSL on Windows)

### Build

```bash
g++ -std=c++17 -o analyzer main.cpp fasta.cpp
```

### Run

```bash
./analyzer
```

### Input Format

The program expects a standard FASTA file:

```
>SequenceID Description
ATGCGATCGATCGATCG
ATCGATCGATCG
```

Each record begins with a `>` header line followed by one or more lines of sequence data.

## Codon Table

The program includes a full RNA codon table mapping all 64 codons to their corresponding amino acids, including the three stop codons (`UAA`, `UAG`, `UGA`).

## Roadmap

- [x] FASTA file parsing
- [ ] Sequence statistics (length, composition, molecular weight)
- [ ] Motif search via Trie
- [ ] DNA to protein translation
- [ ] Smith-Waterman local alignment

## Contributing

This is a personal learning project. Feel free to open an issue or pull request if you spot something worth improving.

## License

MIT License — free to use and modify.
