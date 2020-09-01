Scripts
=======

Some useful scripts to handle fasta files.

- **BreakLines.py**: script designed to break huge sequences in a single line to 100 nucleotides length per line.
    - Usage: ```BreakLines.py input.fasta output_breaklines.fasta```

- **RemoveRedundancy.py**: script designed to remove redundancy from datasets (cluster sequences 100% similar).
    - Usage: ```RemoveRedundancy.py input.fa output_RedundancyRemoved.fa report.txt```

- **TranslatePartial.py**: script designed to translate the CDSs predicted by CodAn using the PARTIAL models.
    - Usage: ```TranslatePartial.py input.fa output.fa```
    
