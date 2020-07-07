


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

# Perform HQ read extraction


inputFile = inputReadsFile
outfile = open(outputFolder+"/hq_reads.fasta","w")
numSeq = 0
for a in range(256,55,-50):
	print("Extracting high confident reads with length equal to %d" %a)
	totBases = 0
	os.system("kmc -k"+str(a)+" "+inputFile+" "+outputFolder+"/kmcOutput "+outputFolder+"/")
	os.system("kmc_dump -ci3 "+outputFolder+"/kmcOutput "+outputFolder+"/kmcDump_output")
	infile = open(outputFolder+"/kmcDump_output")
	while True:
		line = infile.readline().rstrip()
		if not line:
			break
		fields = line.split("\t")
		numSeq+=1
		outfile.write(">Sequence_"+str(numSeq)+"\n"+fields[0]+"\n")
		totBases+=1
	infile.close()
	print("Found %d high confident reads")
	totBases = totBases*a
	if totBases>10*genomeSize:
		break
outfile.close()





reads = outputFolder+"/hq_reads.fasta"



#Perform reference guided de novo assembly
if reads == "":
	print("Something went wrong with the quality filtering step, now exiting")
	exit()

else:
	
	readsSeq = {}
	
	
	print("* Reference guided de novo assembly")

	

	stage_a = open(outputFolder+"/local_assemblies.fasta","w")




	print("* * Assembly on sliding windows started")
	
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
		print("Aligning hq reads with bowtie2")
		os.system(installationDirectory+"/src/conda/bin/bowtie2-build "+outputFolder+"/partReference.fasta "+outputFolder+"/reference "+outputFolder+"/null")
		os.system(installationDirectory+"/src/conda/bin/bowtie2  -p "+numThreads+" -x "+outputFolder+"/reference -f  "+outputFolder+"/hq_reads.fasta -S "+outputFolder+"/alignment.sam")
		print("Converting to bam")
		os.system(installationDirectory+"/src/conda/bin/samtools view -bS -h -F 4 "+outputFolder+"/alignment.sam > "+outputFolder+"/alignment.bam")
		#os.system(installationDirectory+"/src/conda/bin/samtools sort -o "+outputFolder+"/alignment_sorted.bam "+outputFolder+"/alignment.bam")
		print("Extracting aligned reads")
		os.system("bam2fastq -o "+outputFolder+"/alignedReads -f -q "+outputFolder+"/alignment.bam")
		print("Performing local assembly")
		sys.stdin.read(1)
		os.system("spades.py -s "+outputFolder+"/alignment")
		
		for seq_record in SeqIO.parse(outputFolder+"/outputIdba/scaffold.fa","fasta"):
			if len(str(seq_record.seq)) > maxScaffoldLength:
				maxScaffoldLength = len(str(seq_record.seq))
				longestContig = str(seq_record.seq)

		print("* * * Contig size: "+str(maxScaffoldLength))
			
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