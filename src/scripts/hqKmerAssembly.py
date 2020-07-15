import os,sys
from Bio import SeqIO
from Bio import Seq
import argparse
parser = argparse.ArgumentParser(description="Attempt an assembly with highly frequent kmers")
parser.add_argument("-p","--installationDirectory",required=True,help="The installation directory")
parser.add_argument("-r","--reads",required=True,help="fasta file with reads mapped to a genome or a portion")
parser.add_argument("-ref","--reference",required=True,help="fasta file with a reference genome")
parser.add_argument("-t","--threads",required=True,help="Number of threads")
parser.add_argument("-of","--outputFolder",required=True,help="output folder")
args = vars(parser.parse_args())

installationDirectory = args['installationDirectory']
numThreads = args['threads']
reference  = args['reference']
reads = args['reads']
outputFolder = args['outputFolder']

for seq_record in SeqIO.parse(reference,"fasta"):
    refSequence = str(seq_record.seq)
    refLength = len(refSequence)

#Calculating coverag
coverage = 0
for seq_record in SeqIO.parse(reads,"fasta"):
    coverage += len(str(seq_record.seq))

coverage = float(coverage)/float(refLength)/10
print("Calculated coverage %f" %coverage)
assembledSequence = ""
kmerSize = 81
kmerCoverage = 50
fractionKmer = 1
while float(len(assembledSequence))/float(refLength) < 0.8:
    print("Trying kmer size %d / kmer coverage %f" %(kmerSize,kmerCoverage))
    #prepare kmer database 
    print("kmc -k"+str(kmerSize)+" "+reads+" kmerDB "+outputFolder)
    os.system("kmc -k"+str(kmerSize)+" "+reads+" kmerDB "+outputFolder)
    os.system("kmc_dump -ci"+str(int(kmerCoverage))+" "+outputFolder+"/kmerDB "+outputFolder+"/kmerDB_output")
    infile = open(outputFolder+"/kmerDB_output")
    outfile = open(outputFolder+"/kmerDB_output.fastq","w")
    numSeq = 0
    os.system("rm -rf "+outputFolder+"/outputSpades")
    while True:
        line = infile.readline().rstrip()
        if not line:
            break
        fields = line.split("\t")
        numSeq+=1
        outfile.write("@Sequence_"+str(numSeq)+"\n"+fields[0]+"\n+\n")
        for a in range(len(fields[0])):
            outfile.write("G")
        outfile.write("\n")
    outfile.close()
    os.system("spades.py -s "+outputFolder+"/kmerDB_output.fastq --phred-offset 33 --careful -o "+outputFolder+"/outputSpades")

    maxScaffoldLength = 0

    if os.path.isfile(outputFolder+"/outputIdba/scaffold.fa") == True:
        for seq_record in SeqIO.parse(outputFolder+"/outputSpades/scaffolds.fasta","fasta"):
            if len(str(seq_record.seq)) > maxScaffoldLength:
                maxScaffoldLength = len(str(seq_record.seq))
                assembledSequence = str(seq_record.seq)
    fractionKmer +=1
    kmerCoverage = kmerCoverage - 10
    sys.stdin.read(1)
    

