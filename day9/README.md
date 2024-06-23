# DNA Analyzer

This Python script analyzes DNA sequences from FASTA or GeneBank files to either find the longest repeated subsequence or identify the 50bp region with the highest GC content.

## Features

- **Find Longest Repeated Subsequence**: Identifies the longest subsequence that appears at least twice within the DNA sequence.
- **Find Highest GC Content Area**: Identifies the 50bp region with the highest GC content within the DNA sequence.

## Requirements

- Python 3.x
- Biopython

## usage
python analyze.py <path_to_fasta_or_genbank_file> <function>

Examples:
- **Find Longest Repeated Subsequence**: ~ python3 analyze.py path/to/your/file.fasta duplicate
- **Find Highest GC Content Area**: ~ python3 analyze.py path/to/your/file.fasta high_GC_area
