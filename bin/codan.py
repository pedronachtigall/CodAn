#!/usr/bin/env python3

#CodAn predicts CDS and UTR sequences in Full-Length and partial assembled transcripts from transcriptome data
#Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except:
    print('''CodAn was not able to run due to the ERROR below:
    The biopython package is not properly installed or it is not active in your actual environment.
    Please, install biopython properly (check https://biopython.org/ for more details), or active the environment where biopython is installed!''')
    quit()
from optparse import OptionParser
import datetime as dt
import subprocess
import os

#>>>>Translate ORFs
def _TranslateORFs_(outF):
    translated = open(outF+"PEP_sequences.fa","w")
    for record in SeqIO.parse(outF+"ORF_sequences.fasta", "fasta"):
        ID = str(record.id)
        protein = str(record.seq.translate())
        translated.write(">"+ID+"\n"+protein+"\n")
    translated.close()

#>>>>Get Reverse Complement of transcripts for ORF prediction (used in _predict_)
def _reverse_complement_(transcripts, outF): #MINUS
    minus = open(outF+"minus.fa","w")
    for record in SeqIO.parse(transcripts, "fasta"):
        sequence = str(record.seq.reverse_complement())
        minus.write(">"+str(record.id)+"\n"+"\n".join([sequence[n:n+100] for n in range(0,len(sequence),100)])+"\n")
    minus.close()

#>>>>Read prediction result (used in _retrieveORF_)
def _readORF_(ORF):
    start = []
    stop = []
    cds = []
    a = open(ORF, "r")
    for line in a:
        if "\tCDS\t" in line:
            line1 = line.rstrip().split("\t")
            cds.append([line1[0], line1[3], line1[4]])
        if "\tstop_codon\t" in line:
            line1 = line.rstrip().split("\t")
            stop.append([line1[0], line1[4]])
        if "\tstart_codon\t" in line:
            line1 = line.rstrip().split("\t")
            start.append([line1[0], line1[3]])
    a.close()    
    return cds, stop, start

def _readORF_BOTH_(ORF):
    start = {}
    stop = {}
    cds = {}
    a = open(ORF, "r")
    for line in a:
        if "\tCDS\t" in line:
            line1 = line.rstrip().split("\t")
            cds[line1[0]] = [line1[3], line1[4], line1[5]]
        if "\tstop_codon\t" in line:
            line1 = line.rstrip().split("\t")
            stop[line1[0]] = [line1[3], line1[4], line1[5]]
        if "\tstart_codon\t" in line:
            line1 = line.rstrip().split("\t")
            start[line1[0]] = [line1[3], line1[4], line1[5]]
    a.close()    
    return cds, stop, start

#>>>>Retrieve ORF and UTR sequences from prediction results (used in _codan_)
def _retrieveORF_PLUS_(transcripts, outF):
    cds, stop, start = _readORF_(outF+"ORFs.gtf")
    record_dictP = SeqIO.index(transcripts, "fasta")

    wORF = open(outF+"ORF_sequences.fasta", "w")
    for i in cds:
        sequence = record_dictP[i[0]].seq
        orf_seq = str(sequence[int(i[1])-1:int(i[2])+3])
        wORF.write(">"+i[0]+"\n"+"\n".join([orf_seq[n:n+100] for n in range(0,len(orf_seq),100)])+"\n")
    wORF.close()
    
    tUTR = open(outF+"3utr_sequences.fasta", "w")
    for i in stop:
        sequence = record_dictP[i[0]].seq
        stop_seq = sequence[int(i[1]):]
        if stop_seq == "":
            stop_seq = "None"
        tUTR.write(">"+i[0]+"\n"+str(stop_seq)+"\n")
    tUTR.close()
    
    fUTR = open(outF+"5utr_sequences.fasta", "w")
    for i in start:
        sequence = record_dictP[i[0]].seq
        start_seq = sequence[:int(i[1])-1]
        if start_seq == "":
            start_seq = "None"
        fUTR.write(">"+i[0]+"\n"+str(start_seq)+"\n")
    fUTR.close()
    
    print("\tnumber of transcripts -> "+str(len(record_dictP)))
    print("\tnumber of predictions -> "+str(len(cds)))
    
    os.remove(transcripts+".index")

def _retrieveORF_MINUS_(transcripts, outF):
    cds, stop, start = _readORF_(outF+"ORFs.gtf")
    record_dictP = SeqIO.index(outF+"minus.fa", "fasta")

    wORF = open(outF+"ORF_sequences.fasta", "w")
    for i in cds:
        sequence = record_dictP[i[0]].seq
        orf_seq = str(sequence[int(i[1])-1:int(i[2])+3])
        wORF.write(">"+i[0]+"\n"+"\n".join([orf_seq[n:n+100] for n in range(0,len(orf_seq),100)])+"\n")
    wORF.close()
    
    tUTR = open(outF+"3utr_sequences.fasta", "w")
    for i in stop:
        sequence = record_dictP[i[0]].seq
        stop_seq = sequence[int(i[1]):]
        if stop_seq == "":
            stop_seq = "None"
        tUTR.write(">"+i[0]+"\n"+str(stop_seq)+"\n")
    tUTR.close()
    
    fUTR = open(outF+"5utr_sequences.fasta", "w")
    for i in start:
        sequence = record_dictP[i[0]].seq
        start_seq = sequence[:int(i[1])-1]
        if start_seq == "":
            start_seq = "None"
        fUTR.write(">"+i[0]+"\n"+str(start_seq)+"\n")
    fUTR.close()

    print("\tnumber of transcripts -> "+str(len(record_dictP)))
    print("\tnumber of predictions -> "+str(len(cds)))

    os.remove(outF+"minus.fa")
    os.remove(outF+"minus.fa.index")

def _retrieveORF_BOTH_(transcripts, minus, outF):
    cdsP, stopP, startP = _readORF_BOTH_(outF+"ORFs_plus.gtf")
    cdsM, stopM, startM = _readORF_BOTH_(outF+"ORFs_minus.gtf")

    record_dictP = SeqIO.index(transcripts, "fasta")
    record_dictM = SeqIO.index(minus, "fasta")
    
    plus_strand = []
    minus_strand = []
    wORF = open(outF+"ORF_sequences.fasta", "w")
    for k in sorted(cdsP.keys()):
        if k in cdsM.keys():
            plus = cdsP[k]
            minus = cdsM[k]
            # sizeP = int(plus[1])-int(plus[0])
            # sizeM = int(minus[1])-int(minus[0])
            # if sizeP >= sizeM:
            #     plus_strand.append(k)
            #     i = cdsP[k]
            #     sequence = record_dictP[k].seq
            #     orf_seq = sequence[int(i[0])-1:int(i[1])+3]
            #     wORF.write(">"+k+"\n"+str(orf_seq)+"\n")
            # if sizeP < sizeM:
            #     minus_strand.append(k)
            #     i = cdsM[k]
            #     sequence = record_dictM[k].seq
            #     orf_seq = sequence[int(i[0])-1:int(i[1])+3]
            #     wORF.write(">"+k+"\n"+str(orf_seq)+"\n")
            viterbiP = float(plus[2])
            viterbiM = float(minus[2])
            if viterbiP >= viterbiM:
                plus_strand.append(k)
                i = cdsP[k]
                sequence = record_dictP[k].seq
                orf_seq = str(sequence[int(i[0])-1:int(i[1])+3])
                wORF.write(">"+k+"\n"+"\n".join([orf_seq[n:n+100] for n in range(0,len(orf_seq),100)])+"\n")
            if viterbiP < viterbiM:
                minus_strand.append(k)
                i = cdsM[k]
                sequence = record_dictM[k].seq
                orf_seq = str(sequence[int(i[0])-1:int(i[1])+3])
                wORF.write(">"+k+"\n"+"\n".join([orf_seq[n:n+100] for n in range(0,len(orf_seq),100)])+"\n")
        if k not in cdsM.keys():
            plus_strand.append(k)
            i = cdsP[k]
            sequence = record_dictP[k].seq
            orf_seq = str(sequence[int(i[0])-1:int(i[1])+3])
            wORF.write(">"+k+"\n"+"\n".join([orf_seq[n:n+100] for n in range(0,len(orf_seq),100)])+"\n")
    for k in sorted(cdsM.keys()):
        if k not in cdsP.keys():
            minus_strand.append(k)
            i = cdsM[k]
            sequence = record_dictM[k].seq
            orf_seq = str(sequence[int(i[0])-1:int(i[1])+3])
            wORF.write(">"+k+"\n"+"\n".join([orf_seq[n:n+100] for n in range(0,len(orf_seq),100)])+"\n")
    wORF.close()
    
    tUTR = open(outF+"3utr_sequences.fasta", "w")
    for k in sorted(stopP.keys()):
        if k in plus_strand:
            i = stopP[k]
            sequence = record_dictP[k].seq
            stop_seq = sequence[int(i[1]):]
            if stop_seq == "":
                stop_seq = "None"
            tUTR.write(">"+k+"\n"+str(stop_seq)+"\n")
    for k in sorted(stopM.keys()):
        if k in minus_strand:
            i = stopM[k]
            sequence = record_dictM[k].seq
            stop_seq = sequence[int(i[1]):]
            if stop_seq == "":
                stop_seq = "None"
            tUTR.write(">"+k+"\n"+str(stop_seq)+"\n")
    tUTR.close()
    
    fUTR = open(outF+"5utr_sequences.fasta", "w")
    for k in sorted(startP.keys()):
        if k in plus_strand:
            i = startP[k]
            sequence = record_dictP[k].seq
            start_seq = sequence[:int(i[0])-1]
            if start_seq == "":
                start_seq = "None"
            fUTR.write(">"+k+"\n"+str(start_seq)+"\n")
    for k in sorted(startM.keys()):
        if k in minus_strand:
            i = startM[k]
            sequence = record_dictM[k].seq
            start_seq = sequence[:int(i[0])-1]
            if start_seq == "":
                start_seq = "None"
            fUTR.write(">"+k+"\n"+str(start_seq)+"\n")
    fUTR.close()
    
    #adjust and mix gtf of plus and minus
    outG = open(outF+"annotation.gtf","w")
    PLUSgtf = open(outF+"ORFs_plus.gtf","r")
    for line in PLUSgtf:
        if not line.startswith("\n"):
            line1 = line.rstrip().split("\t")
            if line1[0] in plus_strand:
                outG.write(line)
    PLUSgtf.close()
    MINUSgtf = open(outF+"ORFs_minus.gtf","r")
    for line in MINUSgtf:
        if not line.startswith("\n"):
            line1 = line.rstrip().split("\t")
            if line1[0] in minus_strand:
                if "\tstart_codon\t" in line:
                    line1 = line.rstrip().split("\t")
                    seq = record_dictP[line1[0]]
                    st = len(seq)-int(line1[4])
                    end = len(seq)-int(line1[3])
                    line1[3] = str(st+1)
                    line1[4] = str(end+1)
                    line1[6] = "-"
                    outG.write("\t".join(line1)+"\n")
                if "\tCDS\t" in line:
                    line1 = line.rstrip().split("\t")
                    seq = record_dictP[line1[0]]
                    st = len(seq)-int(line1[4])
                    end = len(seq)-int(line1[3])
                    line1[3] = str(st+1)
                    line1[4] = str(end+1)
                    line1[6] = "-"
                    outG.write("\t".join(line1)+"\n")
                if "stop_codon" in line:
                    line1 = line.rstrip().split("\t")
                    seq = record_dictP[line1[0]]
                    st = len(seq)-int(line1[4])
                    end = len(seq)-int(line1[3])
                    line1[3] = str(st+1)
                    line1[4] = str(end+1)
                    line1[6] = "-"
                    outG.write("\t".join(line1)+"\n")
    MINUSgtf.close()
    outG.close()
    
    print("\tnumber of transcripts -> "+str(len(record_dictP)))
    print("\tnumber of predictions -> "+str(len(plus_strand)+len(minus_strand)))
    print("\t\tpredictions at plus strand -> "+str(len(plus_strand)))
    print("\t\tpredictions at minus strand -> "+str(len(minus_strand)))
    
    os.remove(outF+"ORFs_plus.gtf")
    os.remove(outF+"ORFs_minus.gtf")
    os.remove(outF+"minus.fa")
    os.system("rm "+outF+"minus.fa.*")
    os.system("rm "+transcripts+".index*")

#>>>>Predict CDS in transcripts (used in _codan_)
def _predict_PLUS_(transcripts, outF, model, cpu):
    # args = ["predict", "-z", model,
    #         "-c", str(cpu),
    #         "-f", transcripts,
    #         ">", outF+"ORFs.gtf"]
    # subprocess.call(args)
    callsequence = "predict -z "+model+" -c "+str(cpu)+" -f "+transcripts+" > "+outF+"ORFs.gtf"+" 2> /dev/null"
    os.system(callsequence)

def _predict_MINUS_(transcripts, outF, model, cpu):
    
    _reverse_complement_(transcripts, outF)
    
    # args = ["predict", "-z", model,
    #         "-c", str(cpu),
    #         "-f", outF+"minus.fa",
    #         ">", outF+"ORFs.gtf"]
    # subprocess.call(args)
    callsequence = "predict -z "+model+" -c "+str(cpu)+" -f "+outF+"minus.fa"+" > "+outF+"ORFs.gtf"+" 2> /dev/null"
    os.system(callsequence)

def _predict_BOTH_(transcripts, outF, model, cpu):
    # args = ["predict", "-z", model,
    #         "-c", str(cpu),
    #         "-f", transcripts,
    #         ">", outF+"ORFs_plus.gtf"]
    # subprocess.call(args)
    callsequence = "predict -z "+model+" -c "+str(cpu)+" -f "+transcripts+" > "+outF+"ORFs_plus.gtf"+" 2> /dev/null"
    os.system(callsequence)

    _reverse_complement_(transcripts, outF)
    
    # args = ["predict", "-z", model,
    #         "-c", str(cpu),
    #         "-f", outF+"minus.fa",
    #         ">", outF+"ORFs_minus.gtf"]
    # subprocess.call(args)
    callsequence = "predict -z "+model+" -c "+str(cpu)+" -f "+outF+"minus.fa"+ " > " +outF+"ORFs_minus.gtf"+" 2> /dev/null"
    os.system(callsequence)

#>>>>Run CodAn
def _codan_PLUS_(transcripts, outF, model, cpu):
    print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> CDS prediction...")
    _predict_PLUS_(transcripts, outF, model, cpu)
    print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> retrieving sequences...")
    _retrieveORF_PLUS_(transcripts, outF)
    if "_full" in model:
        _TranslateORFs_(outF)

def _codan_MINUS_(transcripts, outF, model, cpu):
    print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> CDS prediction...")
    _predict_MINUS_(transcripts, outF, model, cpu)
    print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> retrieving sequences...")
    _retrieveORF_MINUS_(outF+"minus.fa", outF)
    if "_full" in model:
        _TranslateORFs_(outF)

def _codan_BOTH_(transcripts, outF, model, cpu):
    print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> CDS prediction...")
    _predict_BOTH_(transcripts, outF, model, cpu)
    print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> retrieving sequences...")
    _retrieveORF_BOTH_(transcripts, outF+"minus.fa", outF)
    if "_full" in model:
        _TranslateORFs_(outF)

##>>>>BLAST search
def _blast_search_(outF, blastDB, HSP, threads):
    print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> running blast search...")
    
    # args = ["blastx", "-query", outF+"ORF_sequences.fasta",
    #         "-db", blastDB,
    #         "-outfmt", "6",
    #         "-max_target_seqs", "1",
    #         "-qcov_hsp_perc", str(HSP),
    #         "-out", outF+"blast_results/"+"blast_result.tabular"]
    # subprocess.call(args)
    callsequence = "blastx -query "+outF+"ORF_sequences.fasta -db "+blastDB+" -outfmt 6 -max_target_seqs 1 -qcov_hsp_perc "+str(HSP)+" -out "+outF+"blast_results/blast_result.tabular -num_threads "+str(threads)
    os.system(callsequence)
    
    print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> parsing blast results...")
    
    SQ = {}
    a = open(outF+"blast_results/"+"blast_result.tabular", "r")
    for line in a:
        line1 = line.rstrip().split("\t")
        SQ[line1[0]] = line1[1]        
    a.close()
    
    LOG = open(outF+"blast_results/"+"report.log", "w")
    LOG.write("transcriptID\thit_ID\n")
    
    orf_blast = open(outF+"blast_results/ORF_blast_hit.fasta", "w")
    for record in SeqIO.parse(outF+"ORF_sequences.fasta", "fasta"):
        ID = record.id
        sequence = record.seq
        if ID in SQ.keys():
            LOG.write(str(ID)+"\t"+SQ[ID]+"\n")
            orf_blast.write(">"+str(ID)+"__"+SQ[ID]+"\n"+str(sequence)+"\n")
        else:
            LOG.write(str(ID)+"\tno hit found\n")
            orf_blast.write(">"+str(ID)+"__no_hit_found\n"+str(sequence)+"\n")
    orf_blast.close()
    
    LOG.close()

    tUTR_blast = open(outF+"blast_results/3utr_blast_hit.fasta", "w")
    for record in SeqIO.parse(outF+"3utr_sequences.fasta", "fasta"):
        ID = record.id
        sequence = record.seq
        if ID in SQ.keys():
            tUTR_blast.write(">"+str(ID)+"__"+SQ[ID]+"\n"+str(sequence)+"\n")
        else:
            tUTR_blast.write(">"+str(ID)+"__no_hit_found\n"+str(sequence)+"\n")
    tUTR_blast.close()
    
    fUTR_blast = open(outF+"blast_results/5utr_blast_hit.fasta", "w")
    for record in SeqIO.parse(outF+"5utr_sequences.fasta", "fasta"):
        ID = record.id
        sequence = record.seq
        if ID in SQ.keys():
            fUTR_blast.write(">"+str(ID)+"__"+SQ[ID]+"\n"+str(sequence)+"\n")
        else:
            fUTR_blast.write(">"+str(ID)+"__no_hit_found\n"+str(sequence)+"\n")
    fUTR_blast.close()

##>>>>Options
def __main__():
    parser = OptionParser()
    parser.add_option("-t", "--transcripts", dest="transcripts", help="Mandatory - input transcripts file (FASTA format), /path/to/transcripts.fa", metavar="file", default=None)
    parser.add_option("-m", "--model", dest="model", help="Mandatory - path to model, /path/to/model", metavar="model", default=None)
    parser.add_option("-s", "--strand", dest="strand", help="Optional - strand of sequence to predict genes (plus, minus or both) [default=plus]", metavar="string", default="both")
    parser.add_option("-c", "--cpu", dest="cpu", help="Optional - number of threads to be used [default=1]", metavar="int", default="1")
    parser.add_option("-o", "--output", dest="output_folder", help="Optional - path to output folder, /path/to/output/folder/\nif not declared, it will be created at the transcripts input folder [default=\"CodAn_output\"]", metavar="folder", default=None)
    parser.add_option("-b", "--blastdb", dest="blastDB", help="Optional - path to blastDB of known protein sequences, /path/to/blast/DB/DB_name", metavar="proteinDB", default=None)
    parser.add_option("-H", "--HSP", dest="HSP", help="Optional - used in the \"-qcov_hsp_perc\" option of blastx [default=80]", metavar="int", default=80)
        
    (options, args) = parser.parse_args()
    
    if options.transcripts == None or options.model == None:
        print("""
  _____           _  ___        
 /  __ \         | |/ _ \       
 | /  \/ ___   __| / /_\ \_ __  
 | |    / _ \ / _` |  _  | '_ \ 
 | \__/\ (_) | (_| | | | | | | |
  \____/\___/ \__,_\_| |_/_| |_|
                          
                               
>>>> CodAn v1.1 June 2021 <<<<
****Use -h for help!****

BASIC USAGE (find CDS and UTR sequences):
codan.py -t transcripts.fa -o output_folder -m model_folder

ALTERNATIVE USAGE (predict CDS and UTR sequences and perform BLAST search  in specific DB to annotated predicted genes based on similarity):
codan.py -t transcripts.fa -o output_folder -m model_folder -b blast_DB
        """)
        quit()

    if options.output_folder == None:
        options.output_folder = ""
        folder_l = options.transcripts.split("/")
        for i in range(0, len(folder_l)-1):
            options.output_folder += str(folder_l[i])+"/"
        options.output_folder += "CodAn_output/"
    elif options.output_folder.endswith("/") == False:
        options.output_folder += "/"
    
    if os.path.isdir(options.output_folder) == False:
        os.mkdir(options.output_folder)
    
    if options.model != None and os.path.isdir(options.model) == False:
        print('''CodAn was not able to run due to the ERROR below:
        The model indicated is not a directory.
        CodAn expect to have a directory of the model indicated in the \"-m\" option.
        Please, check the models available and follow the instructions stated at https://github.com/pedronachtigall/CodAn/tree/master/models''')
        quit()
    
    if options.transcripts != None and os.path.isfile(options.transcripts) == False:
        print('''CodAn was not able to run due to the ERROR below:
    The transcripts file indicated is not a valid file.
    Please, indicate a valid fasta file to the \"-t\" option.''')
        quit()

    if options.blastDB != None and os.path.isfile(options.blastDB+".phr") == False:
        print('''CodAn was not able to run due to the ERROR below:
    The blastDB indicated is not a valid protein DB.
    Please, indicate a valid protein DB to the \"-b\" option.
    Check the instructions at https://www.ncbi.nlm.nih.gov/books/NBK279671/ to generate your own protein DB, or donwload a valid protein DB from https://ftp.ncbi.nlm.nih.gov/blast/db/''')
        quit()

    if options.transcripts != None and options.model != None:
        print("""
  _____           _  ___        
 /  __ \         | |/ _ \       
 | /  \/ ___   __| / /_\ \_ __  
 | |    / _ \ / _` |  _  | '_ \ 
 | \__/\ (_) | (_| | | | | | | |
  \____/\___/ \__,_\_| |_/_| |_|
                          
        """)
        print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> starting CodAn (v1.1 June 2021)...")
        print("\ttranscript file -> "+options.transcripts)
        print("\tmodel -> "+options.model)
        print("\tstrand prediction -> "+options.strand)
        print("\tnumber of threads -> "+options.cpu)
        
        if options.strand == "plus":
            _codan_PLUS_(options.transcripts, options.output_folder, options.model, options.cpu)
        if options.strand == "minus":
            _codan_MINUS_(options.transcripts, options.output_folder, options.model, options.cpu)
        if options.strand == "both":
            _codan_BOTH_(options.transcripts, options.output_folder, options.model, options.cpu)

    
        print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> prediction finished!")
        print("\tCDS -> "+options.output_folder+"ORF_sequences.fasta")
        if "_full" in options.model:
            print("\tCDS translated -> "+options.output_folder+"PEP_sequences.fasta")
        print("\t3'UTR -> "+options.output_folder+"3utr_sequences.fasta")
        print("\t5'UTR -> "+options.output_folder+"5utr_sequences.fasta")
        print("\tGTF file -> "+options.output_folder+"annotation.gtf")
    
    if options.blastDB != None:
        if os.path.isdir(options.output_folder+"blast_results/") == False:
            os.mkdir(options.output_folder+"blast_results/")
        _blast_search_(options.output_folder, options.blastDB, options.HSP, options.cpu)
        print("\tBlast search results -> "+options.output_folder+"blast_results/")
        print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> blast analysis finished!")

__main__()

#END
