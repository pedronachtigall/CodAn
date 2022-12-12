![codan_logo](/codan_logo.png)

CodAn
=======
<!---[![GitHub Downloads](https://img.shields.io/github/downloads/pedronachtigall/CodAn/total.svg?style=social&logo=github&label=Download)](https://github.com/pedronachtigall/CodAn/releases) -->
[![Latest GitHub release](https://img.shields.io/github/release/pedronachtigall/CodAn.svg)](https://github.com/pedronachtigall/CodAn/releases/latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3403273.svg)](https://doi.org/10.5281/zenodo.3403273)
[![Published in Briefings in Bioinformatics](https://img.shields.io/badge/published%20in-Briefings%20in%20Bioinformatics-blue)](https://doi.org/10.1093/bib/bbaa045)

**CodAn** (**Cod**ing sequence **An**notator) is a computational tool designed to characterize the CDS and UTR regions on transcripts from any Eukaryote species.

Getting Started
=================

# Installation

Download and decompress the [CodAn.tar.gz](https://github.com/pedronachtigall/CodAn/blob/master/CodAn.tar.gz) file and add the bin directory to your PATH:

```
tar -xf CodAn.tar.gz
export PATH=$PATH:path/to/CodAn/bin/
```

or git clone the CodAn repository and add the bin directory to your PATH:
```
git clone https://github.com/pedronachtigall/CodAn.git
export PATH=$PATH:path/to/CodAn/bin/
```

:warning: If the user is using a macOS, please download the "tops-viterbi_decoding" compiled to macOS [here](https://github.com/pedronachtigall/CodAn/blob/master/for_MacOS_users.zip), decompress the file ```unzip for_MacOS_users.zip```, and copy the "tops-viterbi_decoding" to the bin folder ```cp for_MacOS_users/tops-viterbi_decoding path/to/CodAn/bin/```

# Requirements

- [Python3](https://www.python.org/) and [Biopython](https://biopython.org/wiki/Download)
    - ```apt-get install python3-biopython```
- [Perl](https://www.perl.org/), [Bioperl](https://bioperl.org/) and [MCE](https://metacpan.org/release/MCE) (libmce-perl)
    - ```apt-get install bioperl libmce-perl```
- [NCBI-BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279671/) (v2.9.0 or above)

Ensure that all requirements are working properly.

:warning: **Conda environment installation**

If the user wants to install CodAn and all dependencies using [Conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html), follow the steps below:
- Create the environment:
    - ```conda create -n codan_env python=3.7 biopython perl perl-bioperl perl-mce blast```
- Git clone the CodAn repository and add to your PATH:
    - ```git clone https://github.com/pedronachtigall/CodAn.git```
    - ```export PATH=$PATH:path/to/CodAn/bin/```
- It may be needed to apply "execution permission" to all bin executables:
    - ```chmod 777 path/to/CodAn/bin/*```
- Then, run CodAn as described in the "Usage" section.
- To activate the environment to run CodAn just use the command: ```conda activate codan_env```
- To deactivate the environment just use the command: ```conda deactivate```

:warning: **Conda installation**

[![Install with conda](https://img.shields.io/badge/Install%20with-conda-success)](https://anaconda.org/bioconda/codan)

CodAn can be installed with Conda by using the command: `conda install -c bioconda codan`

The user can also create an environment with the command: `conda create -n codan_env -c bioconda codan`. Then, activate the environment `conda activate codan_env` to run CodAn in Linux and MacOS systems.

- Please, notice that the Conda installation of CodAn does not download the models used in predictions. Download the model specific to your usage [here](https://github.com/pedronachtigall/CodAn/tree/master/models) or using `wget` and decompress the model (by using `unzip` or other tool) to be set in the `-m` parameter.

:warning: **Docker installation**

[![Docker build](https://img.shields.io/badge/Docker-build-blue)](https://hub.docker.com/repository/docker/pedronachtigall/codan)

If the user takes advantage of [Docker](https://docs.docker.com/) in its system, we have a pre-built Dockerfile that allows an easy build and containerization of CodAn. Just follow the steps below:
- Git clone CodAn repository (`git clone https://github.com/pedronachtigall/CodAn.git`) and change to CodAn directory (`cd CodAn`)
- Build the container: `docker build -t codan:v1.0 .` (It may take a few minutes)
- In your working directory (the transcript file should be in there), enter in the container shell: `docker run -v $PWD:/project --rm -it codan:v1.0`
- Just run CodAn by indicating one of the models to the `-m` option: `-m /app/CodAn/models/{VERT|INV|PLANTS|FUNGI}_{full|partial}`
- The command line must be similar to `codan.py -t transcripts.fa -m /app/CodAn/models/MODEL_full/ -o output_folder`

The user may also pull CodAn container direct from the Docker repository following the steps below:
- Pull CodAn container: `docker pull pedronachtigall/codan:latest`
- Run CodAn container: `docker run -v $PWD:/project --rm -it pedronachtigall/codan:latest`
    - Please, notice that you should be in the folder containing your transcripts file
- Use a command line similar to `codan.py -t transcripts.fa -m /app/CodAn/models/MODEL_full/ -o output_folder` to run CodAn and perform its predictions

# Predictive models

The predictive models are available in the subfolder ["models"](https://github.com/pedronachtigall/CodAn/tree/master/models). The folder contains all models designed for Eukaryote species (i.e., Fungi, Plants and Animals [Invertebrates and Vertebrates]). The models were designed to be used in Full-Length or Partial transcripts.

Download the model specific to your necessities, as described at the ["models"](https://github.com/pedronachtigall/CodAn/tree/master/models) folder, decompress the model file (using ```unzip model.zip```), and indicate the decompressed model path in the ```-m``` option.

For example, if you are working with Full-Length transcripts generated from any vertebrate species and will perform the CDS prediction using the Vertebrate Full model.
   - Download the [VERT_full](https://github.com/pedronachtigall/CodAn/blob/master/models/VERT_full.zip) model
   - Decompress the model: ```unzip VERT_full.zip```
   - Indicate the decompressed model to the ```-m``` option: ```-m path/to/VERT_full```

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
codan.py -t transcripts.fa -o output_folder -m model
```

Alternative usage (predict CDS and perform BLAST search in specific DB to annotated predicted genes based on similarity):
```
codan.py -t transcripts.fa -o output_folder -m model -b blast_DB
```
To run this optional step, just indicate a specific protein DB mounted using the ```makeblastdb``` function from the NCBI-BLAST approach.
The user can download the pre-mounted protein DBs, such as swissprot (ftp://ftp.ncbi.nlm.nih.gov/blast/db/).

:warning: :warning: If CodAn not produces any prediction, please check the **Q4** in the FAQ section. :warning: :warning:

# Tutorial
Follow the instructions in the quick [tutorial](https://github.com/pedronachtigall/CodAn/tree/master/tutorial) to learn how to use CodAn and interpret the results.

Reference
=========

If you use or discuss CodAn, please cite:

Nachtigall et al. (2020) CodAn: predictive models for precise identification of coding regions in Eukaryotic transcripts. Briefings in bioinformatics. DOI:[https://doi.org/10.1093/bib/bbaa045](https://doi.org/10.1093/bib/bbaa045)


License
=======

[GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)

Contact
=======
:bug::sos::speech_balloon:

To report bugs, to ask for help and to give any feedback, please contact **Pedro G. Nachtigall**: pedronachtigall@gmail.com

Frequently Asked Questions (FAQ)
================================

**[Q1]** What OS do I need to use CodAn?
- We tested CodAn in Ubuntu 16 and 18. Moreover, we tested CodAn in MacOS Mojave and Catalina, but the user need to download the "tops-viterbi_decoding" compiled to macOS to replace the default "tops-viterbi_decoding" as described in the "Installation" section. In this sense, we believe that CodAn should work on any UNIX OS able to have all dependencies necessary to run CodAn.

**[Q2]** How long does codan need to run an analysis with a set of 200,000 sequences?
- By using 1 thread, the estimated time to analyze 200,000 sequences with CodAn is around 53 minutes. If the user has more threads available for use, which can be set with the option ```-c N``` (where N in the number of threads), the processing time will decrease proportionally as the number of threads being used (e.g., if the user has 6 threads available for the analysis [option ```-c 6```], the processing time of 200,000 sequences will be around 16 minutes). The running time measurement was performed using a personal computer (6-Core i7 with 16Gb memory).

**[Q3]** I am trying to run CodAn on a server and I am not the admin. How can I install all dependencies without using ```apt```?
- You can follow the instructions to use [Conda environments](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html).
- OR, try to install the required perl and python modules following the commands below:
    - Install MCE and Bioperl modules through CPAN (it will install the modules locally):
        - ```perl -MCPAN -e shell```
        ```
        cpan> install MCE
        cpan> install MCE::Mutex
        cpan> install Bio::SeqIO
        cpan> install Bio::DB::Fasta
        cpan> exit
        ```
    - Install Biopython module through pip:
        - ```pip install biopython```

**[Q4]** CodAn is making 0 predictions in the test set. What to do?

- Ensure that all modules are working properly
    - go to the CodAn folder and run each module separetly as follow:
    - tops-viterbi-decoding
    ```
    $tops-viterbi_decoding
    tops-viterbi_decoding: ToPS version "master 00f9ed6"
    Allowed options:
          -h [ --help ]         produce help message
          -m [ --model ] arg    a decodable model
          -F [ --fasta ]        use fasta format
    ```
    - predict
    ```
    $predict
    ERROR: missing fasta file name !
    USAGE: predict [-g <genome> | -t <transcriptome> | -z <local transcriptome predictor> | -s <local genome predictor] -f <fasta file> [-c <number of cpu>]
    ```
    - CodAn
    ```
    $codan.py
       _____           _  ___
      /  __ \         | |/ _ \
      | /  \/ ___   __| / /_\ \_ __
      | |    / _ \ / _` |  _  | '_ \
      | \__/\ (_) | (_| | | | | | | |
       \____/\___/ \__,_\_| |_/_| |_|


    >>>> CodAn v1.0 September 2019 <<<<
    ****Use -h for help!****
    
    BASIC USAGE (find CDS and UTR sequences):
    codan.py -t transcripts.fa -o output_folder -m model_folder
    
    ALTERNATIVE USAGE (predict CDS and UTR sequences and perform BLAST search  in specific DB to annotated predicted genes based on similarity):
    codan.py -t transcripts.fa -o output_folder -m model_folder -b blast_DB
    ```
    - if any different messages print at your terminal, you have problems with the dependencies, try to re-install them.

- Ensure that the headers don't have symbols such as ":" or "|" or " "(space).

- The Bio::DB::Fasta library is responsible for creating the .index. It can't process a fasta file with lines containing more than 65,536 characters. So, if you have any large sequence in one unique line, do the following:
    - download the script [BreakLines.py](https://github.com/pedronachtigall/CodAn/blob/master/scripts/BreakLines.py)
    - run BreakLines script: ```python3 BreakLines.py input.fasta output_breaklines.fasta```
    - use the "output_breaklines.fasta" to run CodAn.

**[Q5]** How can I translate the partial CDSs predicted by the CodAn PARTIAL models?

- You can download the script ```TranslatePartial.py``` [here](https://github.com/pedronachtigall/CodAn/blob/master/scripts/).
    ```
    TranslatePartial.py partialCDS.fa partialCDS_peptide.fa
    ````
