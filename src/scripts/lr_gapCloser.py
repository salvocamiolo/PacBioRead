import os,sys
from Bio import SeqIO
from Bio import Seq
import random as rd
import argparse


parser = argparse.ArgumentParser(description="Join two sequences using long reads and an overlap layaout consensus approach")
parser.add_argument("-p","--installationDirectory",required=True,help="Full path to the folder containing the executables (e.g. conda)")
parser.add_argument("-i","--readsFile",required=True,help="The long read file in fasta format")
#parser.add_argument("-s","--firstBit",required=True,help="The fasta file containing the first sequence to start from in closing the gap")
#parser.add_argument("-e","--secondBit",required=True,help="The fasta file containing the second sequence to start from in closing the gap")
parser.add_argument("-o","--sequenceOutputName",required=True,help="The output file name, the .fasta suffix will be added by the program")
parser.add_argument("-x","--outputFolder",required=True,help="The output folder where to write the output file")
parser.add_argument("-ref","--reference",required=True,help="The reference file used during the de novo assembly")
parser.add_argument("-s","--scaffolds",required=True,help="Scaffolds obtained during the de novo assembly step")


args = vars(parser.parse_args())
installationDirectory = args['installationDirectory']

readsFile = args['readsFile']
#firstBit = args['firstBit']o
#secondBit = args['secondBit']
sequenceOutputName = args['sequenceOutputName']
outputFolder = args['outputFolder']
reference = args['reference']
scaffolds = args['scaffolds']


randomFolderName = "sc140875_"+str(rd.randint(1,1000))
os.system("mkdir "+randomFolderName)
os.chdir(randomFolderName)
os.system("ln -s "+readsFile)
#os.system("ln -s "+firstBit)
#os.system("ln -s "+secondBit)
os.system("ln -s "+reference)
os.system("ln -s "+scaffolds)
readsFile = (readsFile.split("/"))[-1]
reference = (reference.split("/"))[-1]
scaffolds = (scaffolds.split("/"))[-1]
#Load scaffold sequences in memory
scaffoldSequences = {}
for seq_record in SeqIO.parse(scaffolds,"fasta"):
    if not str(seq_record.id) in scaffoldSequences:
        scaffoldSequences[str(seq_record.id)] = str(seq_record.seq)

#Mapping the reads
os.system(installationDirectory+"/src/conda/bin/lastz "+reference+" "+scaffolds+" --format=general > scaffoldMapping.txt")
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

scaffoldName = []
scaffoldPosition = []
scaffoldOrientation = []
for item in scaffoldInfo:
    print(item,scaffoldInfo[item])
    scaffoldName.append(item)
    scaffoldPosition.append(scaffoldInfo[item][1])
    scaffoldOrientation.append(scaffoldInfo[item][2])
scaffoldPositionOrdered = sorted(scaffoldPosition)
print(scaffoldPosition)
print(scaffoldPositionOrdered)
print(scaffoldPosition.index(scaffoldPositionOrdered[0]))
outfile = open(sequenceOutputName,"w")

# ***************************************************
# ********** Starting joining scaffolds *************
# ***************************************************
joinedScaffold = ""
for a in range(len(scaffoldPositionOrdered)-1):
    print("Joining scaffold "+scaffoldName[scaffoldPosition.index(scaffoldPositionOrdered[a])]+" and "+scaffoldName[scaffoldPosition.index(scaffoldPositionOrdered[a+1])])
    
    #Load the first scaffold to join. If the previous pair was joinen, then use the previous joining as first scaffold
    firstBitFile = open("firstBitFile.fasta","w")
    if joinedScaffold == "":
        firstScaffoldToJoin = scaffoldName[scaffoldPosition.index(scaffoldPositionOrdered[a])]
        firstScaffoldOrientation = scaffoldOrientation[scaffoldPosition.index(scaffoldPositionOrdered[a])]
        if firstScaffoldOrientation == "+":
            firstBitFile.write(">firstBit\n"+scaffoldSequences[firstScaffoldToJoin]+"\n")
        else:
            firstBitFile.write(">firstBit\n"+Seq.reverse_complement(scaffoldSequences[firstScaffoldToJoin])+"\n")
        firstBitFile.close()
    else:
        firstBitFile.write(">firstBit\n"+joinedScaffold+"\n")
        firstBitFile.close()

    


    secondScaffoldToJoin = scaffoldName[scaffoldPosition.index(scaffoldPositionOrdered[a+1])]
    secondScaffoldOrientation = scaffoldOrientation[scaffoldPosition.index(scaffoldPositionOrdered[a+1])]

    secondBitFile = open("secondBitFile.fasta","w")
    if secondScaffoldOrientation == "+":
        secondBitFile.write(">secondBit\n"+scaffoldSequences[secondScaffoldToJoin]+"\n")
    else:
        secondBitFile.write(">secondBit\n"+Seq.reverse_complement(scaffoldSequences[secondScaffoldToJoin])+"\n")
    secondBitFile.close()
    







    sequences = {}
    #loading read file in memory

    for seq_record in SeqIO.parse(readsFile,"fasta"):
        if not str(seq_record.id) in sequences:
            sequences[str(seq_record.id)] = str(seq_record.seq)


    #Collect sequences of the two portion to join
    for seq_record in SeqIO.parse("firstBitFile.fasta","fasta"):
        firstBitSeq = str(seq_record.seq)
    for seq_record in SeqIO.parse("secondBitFile.fasta","fasta"):
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
    os.system(installationDirectory+"/src/conda/bin/makeblastdb -dbtype nucl -in "+readsFile+ ">null 2>&1")
    os.system(installationDirectory+"/src/conda/bin/makeblastdb -dbtype nucl -in  endSeq.fasta >null 2>&1")


    #Searching reads mapping the last portion of the first bit
    gapClosed = False
    elongNum = 0
    elongedSequence = firstBitSeq
    while gapClosed == False:
        elongNum +=1 
        os.system(installationDirectory+"/src/conda/bin/blastn -query startSeq.fasta -db "+readsFile+" -outfmt 6 >outputBlast.txt")
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
        
        print("Mi fermo")
        sys.stdin.read(1)
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
        print("Mi fermo 2")
        sys.stdin.read(1)
        line = blastOutputFile.readline().rstrip()
        if not line == "":
            fields = line.split("\t")
            queryStart = int(fields[6])
            queryEnd = int(fields[7])
            subjectStart = int(fields[8])
            subjectEnd = int(fields[9])

            if (queryEnd -queryStart) > 200:
                gapClosed = True
                joinedScaffold = elongedSequence[:queryEnd]+secondBitSeq[subjectEnd:]
                #print(joinedScaffold)
                #print("Finito")
                #sys.stdin.read(1)

    
    #outfile.write(">"+scaffoldName[scaffoldPosition.index(scaffoldPositionOrdered[a])]+"_elonged\n"+elongedSequence+"\n")

#outfile.write(">lastContig\n"+secondBitSeq+"\n")
outfile.write(">finalScaffold\n"+joinedScaffold+"\n")
outfile.close()    
os.system("cp "+sequenceOutputName+" "+outputFolder+"/")
os.chdir("../")
#os.system("rm -rf "+randomFolderName)



















