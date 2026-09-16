# FASTA file Sequence Analyzer

A command-line bioinformatics tool written in C++ that parses FASTA files and runs a full analysis pipeline on biological sequences. Supports DNA, RNA, and protein input with automatic sequence type detection.

Tested on real sequences from NCBI, including Homo sapiens IGHV (JX432019.1) and tumor protein p53 (NM_000546.6).

![Demo](./demo.gif)

## Features

- **FASTA Parsing** — Reads and extracts sequences and headers from `.fasta` files, including multi-record files
- **Automatic Sequence Detection** — Identifies whether a sequence is DNA, RNA, or protein and routes it through the appropriate pipeline
- **DNA → Protein Translation** — Transcribes DNA to RNA, then translates to protein using a full 64-codon table with stop codon handling
- **Sequence Statistics** — Computes length, amino acid composition, and molecular weight in Daltons
- **Motif Search** — Searches for one or more sequence motifs simultaneously using a Trie data structure
- **Local Alignment** — Smith-Waterman local alignment against a query sequence
- **Multithreaded Processing** — Analyzes multiple records concurrently using `std::thread` and `std::mutex`

## Project Structure

```
Protein-Sequence-Analyzer/
├── main.cpp            # Entry point, coordinates program flow
├── fasta.h             # Declarations for all structs and functions
├── fasta.cpp           # Full pipeline implementation
├── test.cpp            # Unit tests for all major functions
├── EXAMPLE_1.fasta     # Sample FASTA file (IGHV, Homo sapiens)
└── TEST_MULTI.fasta    # Multi-record test file (DNA + protein)
```

## Getting Started

### Prerequisites

- g++ with C++17 support or later
- A terminal (Linux, macOS, or WSL on Windows)

### Build

```bash
g++ -std=c++17 -pthread -o analyzer main.cpp fasta.cpp
```

### Run

```bash
# Basic analysis
./analyzer <filename.fasta>

# With Smith-Waterman alignment against a query sequence
./analyzer <filename.fasta> <query_sequence>

# With custom motifs
./analyzer <filename.fasta> <query_sequence> <motif1> <motif2> ...
```

### Examples

```bash
# Analyze a FASTA file
./analyzer EXAMPLE_1.fasta

# Align all records against the IGHV N-terminal motif
./analyzer EXAMPLE_1.fasta QVQL

# Search for custom motifs and align against p53 N-terminal
./analyzer TEST_MULTI.fasta MEEPQ MF KL AG
```

### Build and Run Tests

```bash
g++ -std=c++17 -pthread -o test test.cpp fasta.cpp
./test
```

## Input Format

Standard FASTA format — one or more records, each with a header line starting with `>` followed by sequence data:

```
>SequenceID Description
ATGCGATCGATCGATCG
ATCGATCGATCG
```

Multi-record files are supported. Each record is analyzed concurrently on its own thread.

## Pipeline

```
Input FASTA
    ↓
detectSequenceType()     — DNA / RNA / PROTEIN
    ↓
transcribeDNA()          — DNA only: T → U
    ↓
translateRNA()           — DNA/RNA: codons → single-letter amino acids
    ↓
computeStats()           — length, composition, molecular weight
    ↓
searchMotifs()           — trie-based motif search
    ↓
smithWaterman()          — local alignment against query (optional)
```

## Goals Completed

- [x] FASTA file parsing
- [x] Sequence statistics (length, composition, molecular weight)
- [x] Motif search via Trie
- [x] DNA to protein translation
- [x] Smith-Waterman local alignment
- [x] Automatic sequence type detection (DNA / RNA / Protein)
- [x] Multithreaded record analysis
- [x] Command-line interface with optional query and motifs
- [x] Unit tests

## Contributing

This is a personal learning project. Feel free to open an issue or pull request if you spot something worth improving.

## License

MIT License — free to use and modify.
