import os,sys
from Bio import SeqIO
from Bio import Seq
import random as rd
import argparse

parser = argparse.ArgumentParser(description="Join contiguous contigs together")
parser.add_argument("-p","--installationDirectory",required=True,help="Full path to the folder containing the executables (e.g. conda)")
parser.add_argument("-c","--contigs",required=True,help="Scaffolds obtained during the de novo assembly step")
parser.add_argument("-r","--reference",required=True,help="The reference file used during the de novo assembly")
parser.add_argument("-o","--outputFolder",required=True,help="TThe output folder of the produced joined scaffold")
parser.add_argument("-s","--reads",required=True,help="TThe PacBio reads to use for closing the gap")

args = vars(parser.parse_args())
installationDirectory = args['installationDirectory']
contigs = args['contigs']
reference = args['reference']
outputFolder = args['outputFolder']
reads = args['reads']

#Mapping the reads
os.system(installationDirectory+"/src/conda/bin/lastz "+reference+" "+contigs+" --format=general > scaffoldMapping.txt")
#Getting position and strand for best alignment
infile = open("scaffoldMapping.txt")
infile.readline().rstrip()
scaffoldInfo = {}
while True:
    line = infile.readline().rstrip()
    if not line:
        break
    fields = line.split("\t")
    if not fields[6] in scaffoldInfo:
        scaffoldInfo[fields[6]] = [0,"",""]
    if int(fields[0]) > scaffoldInfo[fields[6]][0]:
        scaffoldInfo[fields[6]] = [int(fields[0]),int(fields[4]),fields[7]]


#GEnerating a new fasta file with the contigs in the correct orientation
contigsSeq = {}
outfile = open(contigs+"_oriented.fasta","w")
for seq_record in SeqIO.parse(contigs,"fasta"):
    if not str(seq_record.id) in contigsSeq:
        if scaffoldInfo[str(seq_record.id)][2]=="+":    
            contigsSeq[str(seq_record.id)] = str(seq_record.seq)
        else:
            contigsSeq[str(seq_record.id)] = Seq.reverse_complement(str(seq_record.seq))


#Get the contigs names in the righ order
infile = open(contigs)
orderedContigs = []
while True:
    title = infile.readline().rstrip()
    if not title:
        break
    infile.readline().rstrip()
    title = title.replace(">","")
    orderedContigs.append(title)




#Join adjacent contigs
elongingSequence = ""
numElongedSequences = 0
outfile = open(outputFolder+"/test.fasta","w")
for a in range(len(orderedContigs)-1):
    print("joining "+orderedContigs[a]+" and "+orderedContigs[a+1])
    startSeqFile = open("startSeq.fasta","w")
    if elongingSequence == "":
        startSeqFile.write(">start\n"+contigsSeq[orderedContigs[a]]+"\n")
        startSeq = contigsSeq[orderedContigs[a]]
        elongingSequence = startSeq
    else:
        startSeqFile.write(">start\n"+elongingSequence+"\n")
        startSeq = elongingSequence
    startSeqFile.close()

    endSeqFile = open("endSeq.fasta","w")
    endSeqFile.write(">end\n"+contigsSeq[orderedContigs[a+1]]+"\n")
    endSeqFile.close()

    #Mapping the reads on the two scaffolds
    print("Mapping reads on %s" %orderedContigs[a])
    os.system(installationDirectory+"/src/conda/bin/minimap2 startSeq.fasta "+reads+" > "+outputFolder+"/minimapBit1")
    #Scanning minimap2 output to search useful reads
    infile = open(outputFolder+"/minimapBit1")
    usefulReads1 = set()
    alnLen1 = {}
    while True:
        line = infile.readline().rstrip()
        if not line:
            break
        fields = line.split("\t")
        subjectEnd = float(fields[8])
        if subjectEnd > 0.9*float(fields[6]):
            usefulReads1.add(fields[0])
            if not fields[0] in alnLen1:
                alnLen1[fields[0]] = int(fields[8]-fields[7])
    infile.close()

        
    print("Mapping reads on %s" %orderedContigs[a+1])
    os.system(installationDirectory+"/src/conda/bin/minimap2  endSeq.fasta "+reads+" > "+outputFolder+"/minimapBit2")
    infile = open(outputFolder+"/minimapBit2")
    usefulReads2 = set()
    alnLen2 = {}
    while True:
        line = infile.readline().rstrip()
        if not line:
            break
        fields = line.split("\t")
        subjectStart = float(fields[7])
        if subjectStart < 5000:
            usefulReads2.add(fields[0])
            if not fields[0] in alnLen2:
                alnLen2[fields[0]] = int(fields[8]-fields[7])
    infile.close()
    print(len(usefulReads1))
    print(len(usefulReads2))
    print(len(usefulReads1 & usefulReads2))

    for item in (usefulReads1 & usefulReads2):
        print(item,alnLen1[item],alnLen2[item])


    print("mi fermo")
    sys.stdin.read(1)
