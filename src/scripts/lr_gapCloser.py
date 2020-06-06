import os,sys
from Bio import SeqIO
from Bio import Seq
import random as rd
import argparse


parser = argparse.ArgumentParser(description="Join two sequences using long reads and an overlap layaout consensus approach")
parser.add_argument("-p","--installationDirectory",required=True,help="Full path to the folder containing the executables (e.g. conda)")
parser.add_argument("-i","--readsFile",required=True,help="The long read file in fasta format")
parser.add_argument("-s","--firstBit",required=True,help="The fasta file containing the first sequence to start from in closing the gap")
parser.add_argument("-e","--secondBit",required=True,help="The fasta file containing the second sequence to start from in closing the gap")
parser.add_argument("-o","--sequenceOutputName",required=True,help="The output file name, the .fasta suffix will be added by the program")
parser.add_argument("-x","--outputFolder",required=True,help="The output folder where to write the output file")
args = vars(parser.parse_args())
installationDirectory = args['installationDirectory']

readFile = args['readsFile']
firstBit = args['firstBit']
secondBit = args['secondBit']
sequenceOutputName = args['sequenceOutputName']
outputFolder = args['outputFolder']



randomFolderName = "sc140875_"+str(rd.randint(1,1000))
os.system("mkdir "+randomFolderName)
os.chdir(randomFolderName)
os.system("ln -s "+readFile)
os.system("ln -s "+firstBit)
os.system("ln -s "+secondBit)
readFile = (readFile.split("/"))[-1]
firstBit = (firstBit.split("/"))[-1]
secondBit = (secondBit.split("/"))[-1]


sequences = {}
#loading read file in memory
if ".fastq" in readsFile or ".fq" in readsFile: 
    for seq_record in SeqIO.parse(readFile,"fastq"):
        if not str(seq_record.id) in sequences:
            sequences[str(seq_record.id)] = str(seq_record.seq)
else:
    for seq_record in SeqIO.parse(readFile,"fasta"):
        if not str(seq_record.id) in sequences:
            sequences[str(seq_record.id)] = str(seq_record.seq)


#Collect sequences of the two portion to join
for seq_record in SeqIO.parse(firstBit,"fasta"):
    firstBitSeq = str(seq_record.seq)
for seq_record in SeqIO.parse(secondBit,"fasta"):
    secondBitSeq = str(seq_record.seq)


#get subsequences to join
startSeq = firstBitSeq[-500:]
startSeqFile = open("startSeq.fasta","w")
startSeqFile.write(">startSeq\n"+startSeq+"\n")
startSeqFile.close()

endSeq = secondBitSeq[:500]
endSeqFile = open("endSeq.fasta","w")
endSeqFile.write(">endSeq\n"+endSeq+"\n")
endSeqFile.close()

#crearing a blast database out of the reads
os.system(installationDirectory+"/src/conda/bin/makeblastdb -dbtype nucl -in "+readFile+ ">null 2>&1")
os.system(installationDirectory+"/src/conda/bin/makeblastdb -dbtype nucl -in  endSeq.fasta >null 2>&1")


#Searching reads mapping the last portion of the first bit
gapClosed = False
elongNum = 0
elongedSequence = firstBitSeq
while gapClosed == False:
    elongNum +=1 
    os.system(installationDirectory+"/src/conda/bin/blastn -query startSeq.fasta -db "+readFile+" -outfmt 6 >outputBlast.txt")
    #Scanning the output file to search best candidates 
    blastOutputFile = open("outputBlast.txt")
    bestWalkingSize = 0
    strandWalking = ""
    walkingSeq = ""
    bestQueryEnd = 0
    bestElongingRead = ""
    while True: #read the first 10 lines only
        line = blastOutputFile.readline().rstrip()
        if not line:
            break
        fields = line.split("\t")
        queryStart = int(fields[6])
        queryEnd = int(fields[7])
        subjectStart = int(fields[8])
        subjectEnd = int(fields[9])

        if (queryEnd - queryStart)> 200:
            if subjectStart < subjectEnd: #alignment plus   
                walkingSize = len(sequences[fields[1]]) - subjectEnd
                if walkingSize > bestWalkingSize:
                    bestWalkingSize = walkingSize
                    strandWalking = "+"
                    walkingSeq = sequences[fields[1]][subjectEnd:]
                    bestQueryEnd = queryEnd
                    bestElongingRead = fields[1]
            else: #The alignment is minus
                walkingSize = subjectEnd
                if walkingSize > bestWalkingSize:
                    bestWalkingSize = walkingSize
                    strandWalking = "-"
                    walkingSeq = Seq.reverse_complement(sequences[fields[1]][:subjectEnd])
                    bestQueryEnd = queryEnd
                    bestElongingRead = fields[1]
    
    #print(elongedSequence[:-500]+"\n\n")
    #print(startSeq[:bestQueryEnd]+"\n\n")
    #print(walkingSeq)
    
    #print(bestWalkingSize,strandWalking,bestElongingRead)
    elongedSequence = elongedSequence[:-500]+startSeq[:bestQueryEnd]+walkingSeq
    #print(elongedSequence)
    #print(elongedSequence)
    print("Sequence length is now %d" %len(elongedSequence))
    startSeq = elongedSequence[-500:]
    startSeqFile = open("startSeq.fasta","w")
    startSeqFile.write(">startSeq elongation number "+str(elongNum)+"\n"+startSeq+"\n")
    startSeqFile.close()
    blastOutputFile.close()

    #Check if the end of the gap has been reached
    elongedSeqFile = open("elongedSeq.fasta","w")
    elongedSeqFile.write(">elonged\n"+elongedSequence+"\n")
    elongedSeqFile.close()
    os.system(installationDirectory+"/src/conda/bin/blastn -query elongedSeq.fasta -db endSeq.fasta -outfmt 6 >outputBlast.txt")
    blastOutputFile = open("outputBlast.txt")
    line = blastOutputFile.readline().rstrip()
    if not line == "":
        fields = line.split("\t")
        queryStart = int(fields[6])
        queryEnd = int(fields[7])
        if (queryEnd -queryStart) > 200:
            gapClosed = True
outfile = open(sequenceOutputName,"w")
outfile.write(">"+sequenceOutputName+"\n"+elongedSequence+"\n")
outfile.close()
os.system("cp "+sequenceOutputName+" "+outputFolder+"/")
os.chdir("../")
os.system("rm -rf "+randomFolderName)



















