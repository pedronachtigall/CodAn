#!/usr/bin/env python3
#script designed to remove redundancy from datasets (cluster sequences 100% similar)
#Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)

import sys
from Bio import SeqIO

def _ParseFastaRR_(fasta):
    final = {}
    for record in SeqIO.parse(fasta,"fasta"):
        ID = str(record.id)
        SEQ = str(record.seq)
        final[ID] = SEQ
    return final

def _ParseFastaInv_(fasta):
    final = {}
    for record in SeqIO.parse(fasta,"fasta"):
        ID = str(record.id)
        SEQ = str(record.seq)
        if SEQ in final.keys():
            final[SEQ].append(ID)
        if SEQ not in final.keys():
            final[SEQ] = [ID]
    return final

def _RemoveRedundancy_(fasta, output, report):
    FASTAinv = _ParseFastaInv_(fasta)
    FASTA = _ParseFastaRR_(fasta)
    R = open(report,"w")
    OUT = open(output,"w")
    countSEQ = 0
    for k in FASTAinv.keys():
        alfa = FASTAinv[k]
        ID = alfa[0]
        SEQ = FASTA[ID]
        seq = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
        OUT.write(">"+ID+"\n"+seq+"\n")
        R.write("//Cluster_"+str(countSEQ)+"\n")
        countID = 0
        for i in alfa:
            R.write(str(countID)+"\t"+i+"\n")
            countID += 1
        countSEQ += 1
    OUT.close()
    R.close()

def _main_():

    if len (sys.argv) != 4:
        print("Script designed to remove redundancy from datasets (cluster sequences 100% similar)")
        print("Basic usage: RemoveRedundancy.py input.fa output.fa report.txt")
        print("\t> input.fa: input file in fasta format")
        print("\t> output.fa: output sequences with break lines in fasta format [100 nts per line]")
        print("\t> report.txt: report file")
        quit()

    fasta = sys.argv[1]
    output = sys.argv[2]
    report = sys.argv[3]
    _RemoveRedundancy_(fasta, output, report)

_main_()

#END