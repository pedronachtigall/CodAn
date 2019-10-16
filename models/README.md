Models
======
This folder contain the models used for CDS prediction:
- **FUNGI_full**: model to predict CDS in Full-Length transcripts of FUNGI species.
- **FUNGI_partial**: model to predict CDS in partial transcripts of FUNGI species.
- **INV_full**: model to predict CDS in Full-Length transcripts of INVERTEBRATE species.
- **INV_partial**: model to predict CDS in partial transcripts of INVERTEBRATE species.
- **PLANTS_full**: model to predict CDS in Full-Length transcripts of PLANT species.
- **PLANTS_partial**: model to predict CDS in partial transcripts of PLANT species.
- **VERT_full**: model to predict CDS in Full-Length transcripts of VERTEBRATE species.
- **VERT_partial**: model to predict CDS in partial transcripts of VERTEBRATE species.

The models are a set of subfolders and files with all parameters specific to each group used in the prediction step of CodAn.


Usage
=====

Download the model specific to your necessities and decompress the file:
```
unzip GROUP_model.zip
```

Indicate the decompressed model to the ```-m``` option of the CodAn tool (```-m path/to/the/GROUP_model```).
