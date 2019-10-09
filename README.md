![codan_logo](/codan_logo.png)

CodAn
=======
<!---[![GitHub Downloads](https://img.shields.io/github/downloads/pedronachtigall/CodAn/total.svg?style=social&logo=github&label=Download)](https://github.com/pedronachtigall/CodAn/releases) -->
[![Latest GitHub release](https://img.shields.io/github/release/pedronachtigall/CodAn.svg)](https://github.com/pedronachtigall/CodAn/releases/latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3403273.svg)](https://doi.org/10.5281/zenodo.3403273)
<!---[![Published in Genome Biology](https://img.shields.io/badge/published%20in-Genome%20Biology-blue.svg)](https://doi.org/10.1101/gr.214270.116) -->

**CodAn** (**Cod**ing sequence **An**notator) is a computational tool designed to characterize the CDS and UTR regions on transcripts from any Eukaryote species.

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

- [Python3](https://www.python.org/) and [Biopython](https://biopython.org/wiki/Download)
    - ```apt-get install python3-biopython```
- [Perl](https://www.perl.org/), [Bioperl](https://bioperl.org/) and [MCE](https://metacpan.org/release/MCE) (libmce-perl)
    - ```apt-get install bioperl libmce-perl```
- [NCBI-BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279671/) (v2.9.0 or above)

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

Basic usage (predict CDS):
```
codan.py -t transcripts.fa -o output_folder -m model_folder
```

Alternative usage (predict CDS and perform BLAST search in specific DB to annotated predicted genes based on similarity):
```
codan.py -t transcripts.fa -o output_folder -m model_folder -b blast_DB
```
To run this optional step, just indicate a specific protein DB mounted using the "makeblastdb" from the NCBI-BLAST approach.
The user can download the pre-mounted protein DBs, such as swissprot (ftp://ftp.ncbi.nlm.nih.gov/blast/db/).

# Predictive models

The predictive models are available in the folder "models". The folder contains all models designed for Eukaryote species (i.e., Fungi, Plants and Animals [Invertebrates and Vertebrates]). The models were designed to be used in Full-Length or partial transcripts. Download the model specific to your necessities, as described at the "models" folder, decompress the model file, and indicate de decompressed model folder path in the "-m" option.

Reference
=========

If you use or discuss CodAn, please cite the preprint:

[Nachtigall et al. CodAn: predictive models for the characterization of mRNA transcripts in Eukaryotes](https://www.biorxiv.org/content/10.1101/794107v1)

License
=======

[GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)

Contact
=======

To report bugs, to ask for help and to give any feedback, please contact **Pedro G. Nachtigall**: pedronachtigall@gmail.com

