#!/usr/bin/env python3
#script to break the sequences into 100 nucleotides per line
#Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)

import sys
from Bio import SeqIO

def _ParseFasta_(fasta):
    final = {}
    for record in SeqIO.parse(fasta, "fasta"):
        ID = str(record.id)
        SEQ = str(record.seq)
        final[ID] = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
    return final

def _GenOutput_(fasta, output):
    FASTA = _ParseFasta_(fasta)
    OUT = open(output,"w")
    for k in FASTA.keys():
        OUT.write(">"+k+"\n"+FASTA[k]+"\n")
    OUT.close()

def _main_():

    if len (sys.argv) != 3:
        print("Basic usage: BreakLines.py input.fa output.fa")
        print("\t> input.fa: input file in fasta format")
        print("\t> output.fa: output sequences with break lines in fasta format [100 nts per line]")
        quit()

    fasta = sys.argv[1]
    output = sys.argv[2]
    _GenOutput_(fasta, output)

_main_()

#END
