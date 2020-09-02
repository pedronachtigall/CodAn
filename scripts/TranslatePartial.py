#!/usr/bin/env python3
#script to translate partial CDSs predicted by CodAn
#Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)

import sys
from Bio import SeqIO
from Bio.Seq import Seq

def _ParseFasta_(fasta):
    final = {}
    for record in SeqIO.parse(fasta, "fasta"):
        ID = str(record.id)
        SEQ = str(record.seq)
        #final[ID] = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
        final[ID] = SEQ
    return final

def _Translate_(cds, frame):
    if frame == 1:
        coding_dna = Seq(cds)
        translated = str(coding_dna.translate())
        return translated
    if frame == 3:
        coding_dna = Seq(cds)
        frame1 = str(coding_dna.translate())
        frame2 = str(coding_dna[1:].translate())
        frame3 = str(coding_dna[2:].translate())
        translated = {1:frame1,2:frame2,3:frame3}
        return translated

def _GenOutput_(fasta, output):
    F = _ParseFasta_(fasta)
    OUT = open(output,"w")
    for k in F.keys():
        SEQ = F[k]
        if SEQ.startswith("ATG"):
            #tranlate with frame 1
            frame1 = _Translate_(SEQ, 1)
            if frame1.count("*") == 1 and frame1[-1] == "*":
                OUT.write(">"+k+" frame1\n"+frame1+"\n")
            if frame1.count("*") <= 1 and not frame1[-1] == "*":
                translated = _Translate_(SEQ, 3)
                for n in translated.keys():
                    if "*" not in translated[n] or (translated[n].count("*") == 1 and translated[n][-1] == "*"):
                        OUT.write(">"+k+" frame"+str(n)+"\n"+translated[n]+"\n")
            if frame1.count("*") > 1:
                translated = _Translate_(SEQ, 3)
                for n in translated.keys():
                    if "*" not in translated[n] or (translated[n].count("*") == 1 and translated[n][-1] == "*"):
                        OUT.write(">"+k+" frame"+str(n)+"\n"+translated[n]+"\n")
        if not SEQ.startswith("ATG"):
            if SEQ[-3:] in ["TAA", "TAG", "TGA"]:
                #translate with 3 frames and check what ends with *
                translated = _Translate_(SEQ, 3)
                for n in translated.keys():
                    if translated[n].endswith("*"):
                        OUT.write(">"+k+" frame"+str(n)+"\n"+translated[n]+"\n")
            if SEQ[-3:] not in ["TAA", "TAG", "TGA"]:
                #translate with 3 frames and keeps all frames with no * in the middle
                translated = _Translate_(SEQ, 3)
                for n in translated.keys():
                    if "*" not in translated[n]:
                        OUT.write(">"+k+" frame"+str(n)+"\n"+translated[n]+"\n")
    OUT.close()

def _main_():

    if len (sys.argv) != 3:
        print("Basic usage: TranslatePartial.py input.fa output.fa")
        print("\t> input.fa: input CDS predicted in fasta format")
        print("\t> output.fa: output translated CDSs in fasta format")
        quit()

    fasta = sys.argv[1]
    output = sys.argv[2]
    _GenOutput_(fasta, output)

_main_()

#END
