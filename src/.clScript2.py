


import os,sys
from Bio import SeqIO
import numpy as np
import time

import argparse

parser = argparse.ArgumentParser(description="Tool to perform reference guided de novo of PacBio reads")
parser.add_argument("-i","--inputReads",required=True, help="The Pacbio reads file in fastq format")
parser.add_argument('-ref',"--reference",required=True,help="Reference file in fasta format")
parser.add_argument("-o","--outputFolder",required=True,help="The output folder name")
parser.add_argument("-q","--quality",required=True, help="Phred quality threshold")
parser.add_argument("-wsize","--windowSize",required=True, help="Window size")
parser.add_argument("-wstep","--windowStep",required=True, help="Window step")
parser.add_argument("-t","--threads",required=True,help="Number of threads")


args = vars(parser.parse_args())
outputFolder = args['outputFolder']
inputReadsFile = args['inputReads']
refFile = args['reference']
threshold = args['quality']
windowSize = int(args['windowSize'])
windowStep = int(args['windowStep'])
numThreads = args['threads']


#Load reference genome in memory

for seq_record in SeqIO.parse(refFile,"fasta"):
	refSeq = str(seq_record.seq)
	genomeSize = len(refSeq)

readsSeq = {}
#Loading high quality reads in memory
numSeq = 0
print("Loading reads in memory")
for seq_record in SeqIO.parse(inputReadsFile,"fastq"):
	numSeq+=1
	if numSeq%50000==0:
		print("%d sequences loaded...." %numSeq)
	if not str(seq_record.id) in readsSeq:
		readsSeq[str(seq_record.id)] = str(seq_record.seq)

# Perform HQ read extraction


inputFile = inputReadsFile

numSeq = 0
hqKmers = []


print("Extracting high confident kmers")
os.system("kmc -k21 "+inputFile+" "+outputFolder+"/kmcOutput "+outputFolder+"/")
os.system("kmc_dump -ci3 "+outputFolder+"/kmcOutput "+outputFolder+"/kmcDump_output")
infile = open(outputFolder+"/kmcDump_output")
while True:
	line = infile.readline().rstrip()
	if not line:
		break
	fields = line.split("\t")
	numSeq+=1
	hqKmers.append(fields[0])

infile.close()
print("Found %d high quality kmers" %(len(hqKmers)))



for a in range(0,len(refSeq),+windowStep):
#for a in range(1):
	endPos = a+windowSize
	if endPos>len(refSeq):
		endPos=len(refSeq)
		windowSize = len(refSeq) - a

	print("* * * Assembling region "+str(a)+"-"+str(endPos))

	partSeq = refSeq[a:endPos]
	tempFasta = open(outputFolder+"/partReference.fasta","w")
	tempFasta.write(">partReference\n"+partSeq+"\n")
	tempFasta.close()

	os.system(installationDirectory+"/src/conda/bin/minimap2 -x map-pb -t "+numThreads+" "+outputFolder+"/partReference.fasta "+inputReadsFile+" > "+outputFolder+"/outputMinimap")

	os.system("awk '($9/$10)>0.70' "+outputFolder+"/outputMinimap | sort -k2rn,2rn >  "+outputFolder+"/outputMinimap_filtered ")

	readsToAssemble = set()
	numAttempt = 0
	maxScaffoldLength = 0

	tfile = open(outputFolder+"/outputMinimap_filtered")
	while True:
		tline = tfile.readline().rstrip()
		if not tline:
			break
		tfields = tline.split("\t")
		readsToAssemble.add(tfields[0])
	tfile.close()

	outfile = open(outputFolder+"/toAssemble.fasta","w")
	numReadsToAssemble = 0
	for item in readsToAssemble:
		if not item == '':
			numReadsToAssemble+=1
			outfile.write(">Sequence_"+str(numReadsToAssemble)+"\n"+readsSeq[item]+"\n")
	outfile.close()
	print("%d reads aligned to the reference " %numReadsToAssemble)

	#Correcting toAssemble with hqKmers
	numCorrReads = 0
	outfile = open(outputFolder+"/toAssemble_corrected.fasta","w")
	for seq_record in SeqIO.parse(outputFolder+"/toAssemble.fasta","fasta"):
		numCorrReads+=1
		if numCorrReads%100==0:
			print("%d reads corrected " %numCorrReads)
		seqID = str(seq_record.id)
		sequence = str(seq_record.seq)
		correctedSequence = ""
		for b in range(0,len(sequence)-21,+21):
			if sequence[b:b+21] in hqKmers:
				correctedSequence+=sequence[b:b+21]
			else:
				correctedSequence+="NNNNNNNNNNNNNNNNNNNNN"







	print("Performing local assembly of all reads with idba")
	os.system(installationDirectory+"/src/conda/bin/art_illumina -i "+outputFolder+"/toAssemble.fasta -l 150 -f 30 -ss HS25 -o "+outputFolder+"/simulatedReads -p -m 500 -s 50")
	toAssembleFile = open(outputFolder+"/allSimulated.fasta","w")
	os.system(installationDirectory+"/src/conda/bin/fq2fa --merge "+outputFolder+"/simulatedReads1.fq "+outputFolder+"/simulatedReads2.fq "+outputFolder+"/allSimulated.fasta")
	os.system("rm -rf "+outputFolder+"/outputIdba/")

	os.system(installationDirectory+"/src/conda/bin/idba_hybrid  --reference "+outputFolder+"/partReference.fasta -r "+outputFolder+"/toAssemble2.fasta --num_threads "+numThreads+" -o "+outputFolder+"/outputIdba > "+outputFolder+"/null 2>&1")
	maxScaffoldLength = 0
	longestContig = ""

	#os.system("rm -rf "+outputFolder+"/sb*")
	#os.system("scaffold_builder_v2.py -q "+outputFolder+"/raven.fasta -r "+outputFolder+"/partReference.fasta -p "+outputFolder+"/sb")
	for seq_record in SeqIO.parse(outputFolder+"/outputIdba/scaffold.fa","fasta"):
		if len(str(seq_record.seq)) > maxScaffoldLength:
			maxScaffoldLength = len(str(seq_record.seq))
			longestContig = str(seq_record.seq)

	print("* * * Contig size: "+str(maxScaffoldLength))

	sys.stdin.read(1)	
	stage_a.write(">Range_"+str(a)+"_"+str(endPos)+"\n"+longestContig+"\n")


stage_a.close()
os.system("rm -rf "+outputFolder+"/outputMinimap* "+outputFolder+"/partReference.fasta " +\
	outputFolder+"/toAssemble.fasta "+outputFolder+"/simulatedReads* "+outputFolder+\
		"/allSimulated.fasta "+outputFolder+"/outputIdba/ " +outputFolder+"/null")

#Joining contigs
print("* * Joining contigs.... ")

os.system(installationDirectory+"/src/conda/bin/python "+installationDirectory+"/src/scripts/contigsJoiner.py -c "+outputFolder+"/local_assemblies.fasta -r "+refFile+" -p "+installationDirectory+" -o "+outputFolder)

#Check final number of scaffold and attempt gap closure if > 1
print("* * Attempting gap closure.... ")

finalScaffols = SeqIO.to_dict(SeqIO.parse(outputFolder+"/scaffolds.fasta","fasta"))

if len(finalScaffols)>1:
	print("\nAttempting gap closure.... ")
	
	os.system(installationDirectory+"/src/conda/bin/python "+installationDirectory+"/src/scripts/oneReadContigsJoiner.py \
		-p "+installationDirectory+ " -c "+outputFolder+"/scaffolds.fasta -r "+refFile+" -x " + \
			outputFolder+" -s "+ inputReadsFile+" -o scaffolds_gapClosed.fasta")
else:
	os.system("mv "+outputFolder+"/scaffolds.fasta "+outputFolder+"/scaffolds_gapClosed.fasta")


#Final alignment and consensus calling
print("* * Calling consensus.... ")
print("* * * Mapping original reads.... ")

os.system(installationDirectory+"/src/conda/bin/minimap2 -a -x map-pb -t "+numThreads+" "+outputFolder+"/scaffolds_gapClosed.fasta "+inputReadsFile+" > "+outputFolder+"/alignment.sam")
print("* * * Converting sam to bam.... ")

os.system(installationDirectory+"/src/conda/bin/samtools view -F 4 -bS -h "+outputFolder+"/alignment.sam > "+outputFolder+"/alignment.bam")
print("* * * Sorting.... ")

os.system(installationDirectory+"/src/conda/bin/samtools sort -o "+outputFolder+"/alignment_sorted.bam "+outputFolder+"/alignment.bam")
print("* * * Indexing.... ")

os.system(installationDirectory+"/src/conda/bin/samtools index "+outputFolder+"/alignment_sorted.bam")
print("* * * Creating pilleup.... ")

os.system(installationDirectory+"/src/conda/bin/samtools mpileup -f "+outputFolder+ \
	"/scaffolds_gapClosed.fasta "+outputFolder +"/alignment_sorted.bam > "+outputFolder+ \
		"/pileup.txt")

print("* * * Calling variants.... ")

os.system(installationDirectory+"/src/conda/bin/varscan mpileup2cns "+outputFolder+"/pileup.txt --variants --output-vcf --min-avg-qual 0 --strand-filter 0 --min-coverage 5   > "+outputFolder+"/output.vcf")
os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/varscanFilter.py -i "+outputFolder+"/output.vcf -o "+outputFolder+"/output_filtered.vcf -1 "+inputReadsFile+" -g 1  -r "+outputFolder+"/scaffolds_gapClosed.fasta -p "+installationDirectory ) 
os.system(installationDirectory+"/src/conda/bin/bgzip -f -c "+outputFolder+"/output_filtered.vcf > "+outputFolder+"/output.vcf_filtered.vcf.gz")
os.system(installationDirectory+"/src/conda/bin/tabix -f "+outputFolder+"/output.vcf_filtered.vcf.gz")
os.system("cat "+outputFolder+"/scaffolds_gapClosed.fasta | "+installationDirectory+"/src/conda/bin/bcftools consensus "+outputFolder+"/output.vcf_filtered.vcf.gz > "+outputFolder+"/finalAssembly.fasta")

print("\n\n De novo assembly finished!")