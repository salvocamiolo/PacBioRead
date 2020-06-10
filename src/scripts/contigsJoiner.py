import os,sys
from Bio import SeqIO
from Bio import Seq
import random as rd
import argparse


parser = argparse.ArgumentParser(description="Join contiguous contigs together")
parser.add_argument("-p","--installationDirectory",required=True,help="Full path to the folder containing the executables (e.g. conda)")
parser.add_argument("-c","--contigs",required=True,help="Scaffolds obtained during the de novo assembly step")
parser.add_argument("-r","--reference",required=True,help="The reference file used during the de novo assembly")

args = vars(parser.parse_args())
installationDirectory = args['installationDirectory']
contigs = args['contigs']
reference = args['reference']





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
    title = title.replace(">","")
    orderedContigs.append(title)


#Join adjacent contigs
for a in range(len(orderedContigs)-1):
    print("joining "+orderedContigs[a]+" and "+orderedContigs[a+1])
    


    






    




