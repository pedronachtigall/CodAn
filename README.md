![codan_logo](/codan_logo.png)

CodAn
=======

A tool to characterize the CDS and UTR regions on any Eukaryote species.

Getting Started
=================

# Installation

Decompress the tar.gz file:

```
tar -xf CodAn.tar.gz
```

Add the bin directory to your PATH:

```
export PATH:$PATH:path/to/CodAn/bin/
```

# Requirements

- Biopython
- Bioperl
- NCBI-BLAST (v2.9.0 or above)

# Usage

```
Usage: codan.py [options]

Options:
  -h, --help            show this help message and exit
  -t file, --transcripts=file
                        Mandatory - input transcripts file (FASTA format),
                        /path/to/transcripts.fa
  -m model, --model=model
                        Mandatory - path to model, /path/to/model
  -s string, --strand=string
                        Optional - strand of sequence to predict genes (plus,
                        minus or both) [default=both]
  -c int, --cpu=int     Optional - number of threads to be used [default=1]
  -o folder, --output=folder
                        Optional - path to output folder,
                        /path/to/output/folder/ if not declared, it will be
                        created at the transcripts input folder
                        [default="CodAn_output"]
  -b proteinDB, --blastdb=proteinDB
                        Optional - path to blastDB of known protein sequences,
                        /path/to/blast/DB/DB_name
  -H int, --HSP=int     Optional - used in the "-qcov_hsp_perc" option of
                        blastx [default=80]

```

Basic usage (find CDS and UTR sequences):
```
codan.py -t transcripts.fa -o output_folder -m model_folder
```

Alternative usage (predict CDS and UTR sequences and perform BLAST search  in specific DB to annotated predicted genes based on similarity):
```
codan.py -t transcripts.fa -o output_folder -m model_folder -b blast_DB
```
To run this optional step, just indicate a specific protein DB mounted using the "makeblastdb" from the NCBI-BLAST approach.
The user can download the pre-mounted protein DBs, such as swissprot, from ftp://ftp.ncbi.nlm.nih.gov/blast/db/ .

# Predictive models

The predictive models are at the folder named "models" and have the models designed for Eukaryote species (i.e., Fungi, Plants and Animals [Invertebrates and Vertebrates]). The models were designed to be used in Full-Length or partial transcripts. Download the model specific to your necessities, as described at the "models" folder, and indicate de decompressed model folder path in the "-m" option.

Reference
=========

If you use CodAn in your analysis, please cite:

Nachtigall et al., under review
