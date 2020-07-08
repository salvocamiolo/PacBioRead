# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'GUI.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog 
from PyQt5.QtWidgets import QMessageBox

import os,sys
from Bio import SeqIO
import numpy as np
import time
import matplotlib.pyplot as plt



class Ui_Form(object):
	def setupUi(self, Form,installationDirectory):
		Form.setObjectName("Form")
		Form.resize(846, 653)
		self.label = QtWidgets.QLabel(Form)
		self.label.setGeometry(QtCore.QRect(10, 10, 101, 16))
		self.label.setObjectName("label")
		self.projectFolderLineEdit = QtWidgets.QLineEdit(Form)
		self.projectFolderLineEdit.setGeometry(QtCore.QRect(10, 30, 361, 21))
		self.projectFolderLineEdit.setObjectName("projectFolderLineEdit")
		self.projectFolderButton = QtWidgets.QPushButton(Form)
		self.projectFolderButton.setGeometry(QtCore.QRect(390, 30, 113, 21))
		self.projectFolderButton.setObjectName("projectFolderButton")
		self.label_2 = QtWidgets.QLabel(Form)
		self.label_2.setGeometry(QtCore.QRect(10, 90, 101, 16))
		self.label_2.setObjectName("label_2")
		self.referenceButton = QtWidgets.QPushButton(Form)
		self.referenceButton.setGeometry(QtCore.QRect(390, 300, 113, 21))
		self.referenceButton.setObjectName("referenceButton")
		self.label_3 = QtWidgets.QLabel(Form)
		self.label_3.setGeometry(QtCore.QRect(10, 280, 101, 16))
		self.label_3.setObjectName("label_3")
		self.referenceLineEdit = QtWidgets.QLineEdit(Form)
		self.referenceLineEdit.setGeometry(QtCore.QRect(10, 300, 361, 21))
		self.referenceLineEdit.setObjectName("referenceLineEdit")
		self.label_13 = QtWidgets.QLabel(Form)
		self.label_13.setGeometry(QtCore.QRect(10, 370, 101, 16))
		self.label_13.setObjectName("label_13")
		self.windowSizeLineEdit = QtWidgets.QLineEdit(Form)
		self.windowSizeLineEdit.setGeometry(QtCore.QRect(10, 390, 101, 21))
		self.windowSizeLineEdit.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.windowSizeLineEdit.setObjectName("windowSizeLineEdit")
		self.windowStepLineEdit = QtWidgets.QLineEdit(Form)
		self.windowStepLineEdit.setGeometry(QtCore.QRect(140, 390, 101, 21))
		self.windowStepLineEdit.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.windowStepLineEdit.setObjectName("windowStepLineEdit")
		self.label_14 = QtWidgets.QLabel(Form)
		self.label_14.setGeometry(QtCore.QRect(140, 370, 101, 16))
		self.label_14.setObjectName("label_14")
		self.runButton = QtWidgets.QPushButton(Form)
		self.runButton.setGeometry(QtCore.QRect(600, 620, 113, 21))
		self.runButton.setObjectName("runButton")
		self.exitButton = QtWidgets.QPushButton(Form)
		self.exitButton.setGeometry(QtCore.QRect(720, 620, 113, 21))
		self.exitButton.setObjectName("exitButton")
		self.label_16 = QtWidgets.QLabel(Form)
		self.label_16.setGeometry(QtCore.QRect(20, 420, 111, 16))
		self.label_16.setObjectName("label_16")
		self.logTextEdit = QtWidgets.QTextEdit(Form)
		self.logTextEdit.setGeometry(QtCore.QRect(10, 440, 821, 161))
		self.logTextEdit.setObjectName("logTextEdit")
		self.numThreadsLineEdit = QtWidgets.QLineEdit(Form)
		self.numThreadsLineEdit.setGeometry(QtCore.QRect(270, 390, 101, 21))
		self.numThreadsLineEdit.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.numThreadsLineEdit.setObjectName("numThreadsLineEdit")
		self.numThreadsLabel = QtWidgets.QLabel(Form)
		self.numThreadsLabel.setGeometry(QtCore.QRect(270, 370, 101, 16))
		self.numThreadsLabel.setObjectName("numThreadsLabel")
		self.label_4 = QtWidgets.QLabel(Form)
		self.label_4.setGeometry(QtCore.QRect(520, 10, 101, 16))
		self.label_4.setObjectName("label_4")
		self.frame = QtWidgets.QFrame(Form)
		self.frame.setGeometry(QtCore.QRect(520, 30, 311, 391))
		self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
		self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
		self.frame.setObjectName("frame")
		self.groupBox = QtWidgets.QGroupBox(self.frame)
		self.groupBox.setGeometry(QtCore.QRect(10, 10, 291, 171))
		self.groupBox.setObjectName("groupBox")
		self.o_numReadsLineEdit = QtWidgets.QLineEdit(self.groupBox)
		self.o_numReadsLineEdit.setEnabled(True)
		self.o_numReadsLineEdit.setGeometry(QtCore.QRect(10, 50, 113, 23))
		self.o_numReadsLineEdit.setReadOnly(True)
		self.o_numReadsLineEdit.setObjectName("o_numReadsLineEdit")
		self.label_5 = QtWidgets.QLabel(self.groupBox)
		self.label_5.setGeometry(QtCore.QRect(10, 30, 101, 16))
		self.label_5.setObjectName("label_5")
		self.o_avQualLineEdit = QtWidgets.QLineEdit(self.groupBox)
		self.o_avQualLineEdit.setEnabled(True)
		self.o_avQualLineEdit.setGeometry(QtCore.QRect(160, 50, 113, 23))
		self.o_avQualLineEdit.setReadOnly(True)
		self.o_avQualLineEdit.setObjectName("o_avQualLineEdit")
		self.label_6 = QtWidgets.QLabel(self.groupBox)
		self.label_6.setGeometry(QtCore.QRect(160, 30, 101, 16))
		self.label_6.setObjectName("label_6")
		self.o_estCoverageLineEdit = QtWidgets.QLineEdit(self.groupBox)
		self.o_estCoverageLineEdit.setEnabled(True)
		self.o_estCoverageLineEdit.setGeometry(QtCore.QRect(10, 120, 113, 23))
		self.o_estCoverageLineEdit.setReadOnly(True)
		self.o_estCoverageLineEdit.setObjectName("o_estCoverageLineEdit")
		self.label_7 = QtWidgets.QLabel(self.groupBox)
		self.label_7.setGeometry(QtCore.QRect(10, 100, 101, 16))
		self.label_7.setObjectName("label_7")
		self.groupBox_2 = QtWidgets.QGroupBox(self.frame)
		self.groupBox_2.setGeometry(QtCore.QRect(10, 210, 291, 171))
		self.groupBox_2.setObjectName("groupBox_2")
		self.r_numReadsLineEdit_3 = QtWidgets.QLineEdit(self.groupBox_2)
		self.r_numReadsLineEdit_3.setEnabled(True)
		self.r_numReadsLineEdit_3.setGeometry(QtCore.QRect(10, 50, 113, 23))
		self.r_numReadsLineEdit_3.setReadOnly(True)
		self.r_numReadsLineEdit_3.setObjectName("r_numReadsLineEdit_3")
		self.label_17 = QtWidgets.QLabel(self.groupBox_2)
		self.label_17.setGeometry(QtCore.QRect(10, 30, 101, 16))
		self.label_17.setObjectName("label_17")
		self.r_avQualLineEdit = QtWidgets.QLineEdit(self.groupBox_2)
		self.r_avQualLineEdit.setEnabled(True)
		self.r_avQualLineEdit.setGeometry(QtCore.QRect(160, 50, 113, 23))
		self.r_avQualLineEdit.setReadOnly(True)
		self.r_avQualLineEdit.setObjectName("r_avQualLineEdit")
		self.label_18 = QtWidgets.QLabel(self.groupBox_2)
		self.label_18.setGeometry(QtCore.QRect(160, 30, 101, 16))
		self.label_18.setObjectName("label_18")
		self.r_estCoverageLineEdit = QtWidgets.QLineEdit(self.groupBox_2)
		self.r_estCoverageLineEdit.setEnabled(True)
		self.r_estCoverageLineEdit.setGeometry(QtCore.QRect(10, 120, 113, 23))
		self.r_estCoverageLineEdit.setReadOnly(True)
		self.r_estCoverageLineEdit.setObjectName("r_estCoverageLineEdit")
		self.label_19 = QtWidgets.QLabel(self.groupBox_2)
		self.label_19.setGeometry(QtCore.QRect(10, 100, 101, 16))
		self.label_19.setObjectName("label_19")
		self.frame_2 = QtWidgets.QFrame(Form)
		self.frame_2.setGeometry(QtCore.QRect(10, 110, 501, 131))
		self.frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
		self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
		self.frame_2.setObjectName("frame_2")
		self.readsFileLineEdit = QtWidgets.QLineEdit(self.frame_2)
		self.readsFileLineEdit.setGeometry(QtCore.QRect(10, 10, 351, 21))
		self.readsFileLineEdit.setObjectName("readsFileLineEdit")
		self.readsFileButton = QtWidgets.QPushButton(self.frame_2)
		self.readsFileButton.setGeometry(QtCore.QRect(380, 10, 113, 21))
		self.readsFileButton.setObjectName("readsFileButton")
		self.qualityStatsButton = QtWidgets.QPushButton(self.frame_2)
		self.qualityStatsButton.setGeometry(QtCore.QRect(10, 50, 241, 21))
		self.qualityStatsButton.setObjectName("qualityStatsButton")
		self.filterButton = QtWidgets.QPushButton(self.frame_2)
		self.filterButton.setGeometry(QtCore.QRect(140, 100, 113, 21))
		self.filterButton.setObjectName("filterButton")
		self.qualityLineEdit = QtWidgets.QLineEdit(self.frame_2)
		self.qualityLineEdit.setGeometry(QtCore.QRect(10, 100, 101, 21))
		self.qualityLineEdit.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.qualityLineEdit.setObjectName("qualityLineEdit")
		self.label_9 = QtWidgets.QLabel(self.frame_2)
		self.label_9.setGeometry(QtCore.QRect(10, 80, 101, 16))
		self.label_9.setObjectName("label_9")

		self.retranslateUi(Form)
		QtCore.QMetaObject.connectSlotsByName(Form)


		self.projectFolderButton.clicked.connect(self.selectProjectFolder)
		self.readsFileButton.clicked.connect(self.selectReadFile)
		self.referenceButton.clicked.connect(self.selectReferenceFile)
		self.runButton.clicked.connect(self.runAssembly)
		self.exitButton.clicked.connect(self.exitProgram)
		self.qualityStatsButton.clicked.connect(self.calculateStatistics)
		self.filterButton.clicked.connect(self.filterReads)
		

	def filterReads(self):
		self.logTextEdit.append("* * Extracting high quality read portions....\n")
		self.logTextEdit.repaint()

		inputFile = self.readsFileLineEdit.text()
		threshold = self.qualityLineEdit.text()
		minLen = 150

		numSeq = 0
		totSequences = 0
		
		totNumBases = 0
		outfile = open(outputFolder+"/hq_reads.fasta","w")
		for seq_record in SeqIO.parse(inputFile,"fastq"):
			numSeq+=1
			totSequences+=1
			if totSequences%3000 == 0:
				self.logTextEdit.append("* * * "+str(totSequences)+" reads analyzed....")
				self.logTextEdit.repaint()
			sequence = str(seq_record.seq)
			totNumBases+=len(sequence)

			quality = seq_record.letter_annotations["phred_quality"]
			maskedSeq = ""
			qualityValues = []
			for a in range(len(quality)):
				qualityValues.append(float(quality[a]))
				if quality[a]>int(threshold):
					maskedSeq+=sequence[a]
				else:
					if len(maskedSeq)>int(minLen):
						outfile.write(">HQSequence_"+str(numSeq)+"\n"+maskedSeq+"\n")
					numSeq+=1
					maskedSeq=""
			if len(qualityValues)==len(quality):
				outfile.write(">HQSequence_"+str(numSeq)+"\n"+maskedSeq+"\n")
				numSeq+=1
				maskedSeq=""

		outfile.close()

		

	
	def calculateStatistics(self):
		inputFile = self.readsFileLineEdit.text()
		threshold = self.qualityLineEdit.text()
		minLen = 150

		self.logTextEdit.append("Read statistics calculation started")
		self.logTextEdit.repaint()
		outputFolder = self.projectFolderLineEdit.text()

		#Load reference genome in memory
		self.logTextEdit.append("Loading reference....")
		self.logTextEdit.repaint()
		refFile = self.referenceLineEdit.text()
		for seq_record in SeqIO.parse(refFile,"fasta"):
			refSeq = str(seq_record.seq)
		
		#Load the reads in memory
		self.logTextEdit.append("Calculating original read stats....")
		self.logTextEdit.repaint()
		
		numSeq = 0
		totSequences = 0
		
		totNumBases = 0
		inputSequences = {}
		outfile = open(outputFolder+"/masked.fasta","w")
		qualityValues = []
		for seq_record in SeqIO.parse(inputFile,"fastq"):
			if not str(seq_record.id) in inputSequences:
				inputSequences[str(seq_record.id)] = seq_record
			numSeq+=1
			totSequences+=1
			
			if totSequences%5000 == 0:
				self.logTextEdit.append("* "+str(totSequences)+" reads analyzed....")
				self.logTextEdit.repaint()
			sequence = str(seq_record.seq)
			totNumBases+=len(sequence)

			quality = seq_record.letter_annotations["phred_quality"]
			if len(qualityValues)<1000000:
				qualityValues+=quality
			
				
		outfile.close()

		self.logTextEdit.append("Calculating stats....")
		self.logTextEdit.repaint()
		self.o_numReadsLineEdit.setText(str(totSequences))
		self.o_avQualLineEdit.setText(str(int(np.mean(qualityValues))))
		self.o_estCoverageLineEdit.setText(str( int(float(totNumBases) / float(len(refSeq))  )  )+" X")
		fig = plt.figure(figsize=(10,10))
		plt.hist(qualityValues,bins=200,label="All reads",color="red",alpha=0.5)
		
	


		self.logTextEdit.append("Calculating reference homologous read stats....")
		self.logTextEdit.append("* Aligning reads....")
		self.logTextEdit.repaint()

		os.system(installationDirectory+"/src/conda/bin/minimap2 -x map-pb -t "+self.numThreadsLineEdit.text()+" "+refFile+" "+inputFile+" > "+outputFolder+"/outputMinimap")
		outfile.close()

		outfile = open(outputFolder+"/refSpecificReads.fastq","w")
		infile = open(outputFolder+"/outputMinimap")
		while True:
			line = infile.readline().rstrip()
			if not line:
				break
			fields = line.split("\t")
			if (float(fields[3])-float(fields[2]))/float(fields[1]) >0.7: #To change from the GUI
				SeqIO.write(inputSequences[fields[0]],outfile,"fastq")
		outfile.close()

		inputFile = outputFolder+"/refSpecificReads.fastq"
		for seq_record in SeqIO.parse(inputFile,"fastq"):

			numSeq+=1
			totSequences+=1
			
			if totSequences%5000 == 0:
				self.logTextEdit.append("* "+str(totSequences)+" reads analyzed....")
				self.logTextEdit.repaint()
			sequence = str(seq_record.seq)
			totNumBases+=len(sequence)

			quality = seq_record.letter_annotations["phred_quality"]
			if len(qualityValues)<1000000:
				qualityValues+=quality
			
				
		outfile.close()

		self.logTextEdit.append("Calculating stats....")
		self.logTextEdit.repaint()
		self.r_numReadsLineEdit_3.setText(str(totSequences))
		self.r_avQualLineEdit.setText(str(int(np.mean(qualityValues))))
		self.r_estCoverageLineEdit.setText(str( int(float(totNumBases) / float(len(refSeq))  )  )+" X")
		plt.hist(qualityValues,bins=200,label="Reference mapping reads",color="blue",alpha=0.5)
		plt.xlabel("Quality phred score")
		plt.legend()
		fig.savefig(outputFolder+"/"+self.readsFileLineEdit.text()+"_qualityDist.png",dpi=300)

		self.logTextEdit.append("\n\nJob finished!")
		self.logTextEdit.repaint()




	def selectProjectFolder(self):
		folderName = QFileDialog.getExistingDirectory(None, "Select project folder","./")
		self.projectFolderLineEdit.setText(folderName)

		if os.path.isdir(folderName) == False:
			os.system("mkdir "+folderName)


	def selectReadFile(self):
		filename, __ = QFileDialog.getOpenFileName(None,"Select read file","./")
		self.readsFileLineEdit.setText(filename)


	def selectReferenceFile(self):
		filename, __ = QFileDialog.getOpenFileName(None,"Select reference file","./")
		self.referenceLineEdit.setText(filename)

	def exitProgram(self):
		exit(1)

	# ****************************************************
	# ***************** Main algorithm *******************
	# ****************************************************
	
	#Check required fields
	def runAssembly(self):

		self.logTextEdit.append("Reference guided de novo assembly started\n")
		self.logTextEdit.repaint()

		if os.path.isdir(self.projectFolderLineEdit.text()) == "":
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("A valid project folder was not selected")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("You should select a folder where all the intermediary file will be stored.\n ")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return

		if os.path.isfile(self.readsFileLineEdit.text()) == False:
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("A valid read file was not selected")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("You should select a fastq file with the PacBio reads.\n ")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return

		if os.path.isfile(self.referenceLineEdit.text()) == False:
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("A valid reference file was not selected")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("You should select a fasta formatted file reference genome.\n ")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return



		#Load reference genome in memory
		refFile = self.referenceLineEdit.text()
		for seq_record in SeqIO.parse(refFile,"fasta"):
			refSeq = str(seq_record.seq)

		
		reads = self.readsFileLineEdit.text()
		if ".fastq" in reads or ".fq" in reads:
			self.logTextEdit.append("* * Converting input file from fastq to fasta....\n")
			self.logTextEdit.repaint()
			outputFolder = self.projectFolderLineEdit.text()
			os.system(installationDirectory+"/src/conda/bin/fq2fa "+self.readsFileLineEdit.text()+" "+outputFolder+'/originalReads.fasta')

		reads = outputFolder+'/originalReads.fasta'
		
		
		#Perform reference guided de novo assembly
		if reads == "":
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("Something did not work with the quality filtering step!")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("The quality filtering step did not produce the expected hq_reads.fasta file in the output folder\n ")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return

		else:
			
			readsSeq = {}
			windowSize = int(self.windowSizeLineEdit.text())
			windowStep = int(self.windowStepLineEdit.text())

			self.logTextEdit.append("* Reference guided de novo assembly")
			self.logTextEdit.append("* * Loading reads in memory")
			self.logTextEdit.repaint()

			#Loading high quality reads in memory
			for seq_record in SeqIO.parse(reads,"fasta"):
				if not str(seq_record.id) in readsSeq:
					readsSeq[str(seq_record.id)] = str(seq_record.seq)

			

			stage_a = open(outputFolder+"/local_assemblies.fasta","w")

			"""self.logTextEdit.append("* * Converting high quality reads to fastq format....")
			self.logTextEdit.repaint()
			with open(reads, "r") as fasta, open(outputFolder+"/hq_reads.fastq", "w") as fastq:
				for record in SeqIO.parse(fasta, "fasta"):
					record.letter_annotations["phred_quality"] = [40] * len(record)
					SeqIO.write(record, fastq, "fastq")"""


			self.logTextEdit.append("* * Assembly on sliding windows started")
			self.logTextEdit.repaint()
			for a in range(0,len(refSeq),+windowStep):
			#for a in range(1):
				endPos = a+windowSize
				if endPos>len(refSeq):
					endPos=len(refSeq)
					windowSize = len(refSeq) - a

				self.logTextEdit.append("* * * Assembling region "+str(a)+"-"+str(endPos))
				self.logTextEdit.repaint()
				partSeq = refSeq[a:endPos]
				tempFasta = open(outputFolder+"/partReference.fasta","w")
				tempFasta.write(">partReference\n"+partSeq+"\n")
				tempFasta.close()

				os.system(installationDirectory+"/src/conda/bin/minimap2 -x map-pb -t "+self.numThreadsLineEdit.text()+" "+outputFolder+"/partReference.fasta "+reads+" > "+outputFolder+"/outputMinimap")

				os.system("awk '($11/$2)>0.70' "+outputFolder+"/outputMinimap | sort -k2rn,2rn >  "+outputFolder+"/outputMinimap_filtered ")

				readsToAssemble = set()
				numAttempt = 0
				maxScaffoldLength = 0

				
				while float(maxScaffoldLength) < float(windowSize)*0.9:
					numAttempt +=1
					if numAttempt == 5:
						break


					for b in range(0,windowSize-500,+150):
						tfile = open(outputFolder+"/outputMinimap_filtered")						

						collectedReads = 0
						while True:
							tline = tfile.readline().rstrip()
							if not tline:
								break
							tfields = tline.split("\t")
							if int(tfields[7]) >b and int(tfields[7]) <(b+150):
								readsToAssemble.add(tfields[0])
								print(tfields[0])
								collectedReads+=1
								if collectedReads == numAttempt:
									break
					tfile.close()

					outfile = open(outputFolder+"/toAssemble.fasta","w")
					numReadsToAssemble = 0
					for item in readsToAssemble:
						if not item == '':
							numReadsToAssemble+=1
							outfile.write(">Sequence_"+str(numReadsToAssemble)+"\n"+readsSeq[item]+"\n")
					outfile.close()

					print("Assembling %d reads" %numReadsToAssemble)
					self.logTextEdit.append("* * * Using "+str(numReadsToAssemble)+" reads....")
					self.logTextEdit.repaint()
					os.system(installationDirectory+"/src/conda/bin/art_illumina -i "+outputFolder+"/toAssemble.fasta -l 150 -f 30 -ss HS25 -o "+outputFolder+"/simulatedReads -p -m 500 -s 50")
					toAssembleFile = open(outputFolder+"/allSimulated.fasta","w")
					os.system(installationDirectory+"/src/conda/bin/fq2fa --merge "+outputFolder+"/simulatedReads1.fq "+outputFolder+"/simulatedReads2.fq "+outputFolder+"/allSimulated.fasta")
					os.system("rm -rf "+outputFolder+"/outputIdba/")
					os.system(installationDirectory+"/src/conda/bin/idba_hybrid  --reference "+outputFolder+"/partReference.fasta -r "+outputFolder+"/allSimulated.fasta --num_threads "+self.numThreadsLineEdit.text()+" -o "+outputFolder+"/outputIdba > "+outputFolder+"/null 2>&1")
					maxScaffoldLength = 0
					longestContig = ""

					for seq_record in SeqIO.parse(outputFolder+"/outputIdba/scaffold.fa","fasta"):
						if len(str(seq_record.seq)) > maxScaffoldLength:
							maxScaffoldLength = len(str(seq_record.seq))
							longestContig = str(seq_record.seq)

					self.logTextEdit.append("* * * Contig size: "+str(maxScaffoldLength))
					self.logTextEdit.repaint()
				stage_a.write(">Range_"+str(a)+"_"+str(endPos)+"\n"+longestContig+"\n")


			stage_a.close()
			os.system("rm -rf "+outputFolder+"/outputMinimap* "+outputFolder+"/partReference.fasta " +\
				outputFolder+"/toAssemble.fasta "+outputFolder+"/simulatedReads* "+outputFolder+\
					"/allSimulated.fasta "+outputFolder+"/outputIdba/ " +outputFolder+"/null")

			#Joining contigs
			self.logTextEdit.append("* * Joining contigs.... ")
			self.logTextEdit.repaint()
			os.system(installationDirectory+"/src/conda/bin/python "+installationDirectory+"/src/scripts/contigsJoiner.py -c "+outputFolder+"/local_assemblies.fasta -r "+refFile+" -p "+installationDirectory+" -o "+outputFolder)

			#Check final number of scaffold and attempt gap closure if > 1
			self.logTextEdit.append("* * Attempting gap closure.... ")
			self.logTextEdit.repaint()
			finalScaffols = SeqIO.to_dict(SeqIO.parse(outputFolder+"/scaffolds.fasta","fasta"))

			if len(finalScaffols)>1:
				self.logTextEdit.append("\nAttempting gap closure.... ")
				self.logTextEdit.repaint()
				os.system(installationDirectory+"/src/conda/bin/python "+installationDirectory+"/src/scripts/oneReadContigsJoiner.py \
					-p "+installationDirectory+ " -c "+outputFolder+"/scaffolds.fasta -r "+refFile+" -x " + \
						outputFolder+" -s "+ self.readsFileLineEdit.text()+" -o scaffolds_gapClosed.fasta")
			else:
				os.system("mv "+outputFolder+"/scaffolds.fasta "+outputFolder+"/scaffolds_gapClosed.fasta")
			
			
			#Final alignment and consensus calling
			self.logTextEdit.append("* * Calling consensus.... ")
			self.logTextEdit.append("* * * Mapping original reads.... ")
			self.logTextEdit.repaint()
			os.system(installationDirectory+"/src/conda/bin/minimap2 -a -x map-pb -t "+self.numThreadsLineEdit.text()+" "+outputFolder+"/scaffolds_gapClosed.fasta "+self.readsFileLineEdit.text()+" > "+outputFolder+"/alignment.sam")
			self.logTextEdit.append("* * * Converting sam to bam.... ")
			self.logTextEdit.repaint()
			os.system(installationDirectory+"/src/conda/bin/samtools view -F 4 -bS -h "+outputFolder+"/alignment.sam > "+outputFolder+"/alignment.bam")
			self.logTextEdit.append("* * * Sorting.... ")
			self.logTextEdit.repaint()
			os.system(installationDirectory+"/src/conda/bin/samtools sort -o "+outputFolder+"/alignment_sorted.bam "+outputFolder+"/alignment.bam")
			self.logTextEdit.append("* * * Indexing.... ")
			self.logTextEdit.repaint()
			os.system(installationDirectory+"/src/conda/bin/samtools index "+outputFolder+"/alignment_sorted.bam")
			self.logTextEdit.append("* * * Creating pilleup.... ")
			self.logTextEdit.repaint()
			os.system(installationDirectory+"/src/conda/bin/samtools mpileup -f "+outputFolder+ \
				"/scaffolds_gapClosed.fasta "+outputFolder +"/alignment_sorted.bam > "+outputFolder+ \
					"/pileup.txt")

			self.logTextEdit.append("* * * Calling variants.... ")
			self.logTextEdit.repaint()
			os.system(installationDirectory+"/src/conda/bin/varscan mpileup2cns "+outputFolder+"/pileup.txt --variants --output-vcf --min-avg-qual 0 --strand-filter 0 --min-coverage 5   > "+outputFolder+"/output.vcf")
			os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/varscanFilter.py -i "+outputFolder+"/output.vcf -o "+outputFolder+"/output_filtered.vcf -1 "+self.readsFileLineEdit.text()+" -g 1  -r "+outputFolder+"/scaffolds_gapClosed.fasta -p "+installationDirectory ) 
			os.system(installationDirectory+"/src/conda/bin/bgzip -f -c "+outputFolder+"/output_filtered.vcf > "+outputFolder+"/output.vcf_filtered.vcf.gz")
			os.system(installationDirectory+"/src/conda/bin/tabix -f "+outputFolder+"/output.vcf_filtered.vcf.gz")
			os.system("cat "+outputFolder+"/scaffolds_gapClosed.fasta | "+installationDirectory+"/src/conda/bin/bcftools consensus "+outputFolder+"/output.vcf_filtered.vcf.gz > "+outputFolder+"/finalAssembly.fasta")

			self.logTextEdit.append("\n\n De novo assembly finished!")
			self.logTextEdit.repaint()




	def retranslateUi(self, Form):
		_translate = QtCore.QCoreApplication.translate
		Form.setWindowTitle(_translate("Form", "Form"))
		self.label.setText(_translate("Form", "Project folder"))
		self.projectFolderButton.setText(_translate("Form", "Select"))
		self.label_2.setText(_translate("Form", "Reads file"))
		self.referenceButton.setText(_translate("Form", "Select"))
		self.label_3.setText(_translate("Form", "Reference file"))
		self.label_13.setText(_translate("Form", "Window size"))
		self.windowSizeLineEdit.setText(_translate("Form", "20000"))
		self.windowStepLineEdit.setText(_translate("Form", "10000"))
		self.label_14.setText(_translate("Form", "Window step"))
		self.runButton.setText(_translate("Form", "Run"))
		self.exitButton.setText(_translate("Form", "Exit"))
		self.label_16.setText(_translate("Form", "Log "))
		self.numThreadsLineEdit.setText(_translate("Form", "8"))
		self.numThreadsLabel.setText(_translate("Form", "Num. threads"))
		self.label_4.setText(_translate("Form", "Read statistics"))
		self.groupBox.setTitle(_translate("Form", "Original"))
		self.o_numReadsLineEdit.setText(_translate("Form", "--"))
		self.label_5.setText(_translate("Form", "Read number"))
		self.o_avQualLineEdit.setText(_translate("Form", "--"))
		self.label_6.setText(_translate("Form", "Average  quality"))
		self.o_estCoverageLineEdit.setText(_translate("Form", "--"))
		self.label_7.setText(_translate("Form", "Coverage"))
		self.groupBox_2.setTitle(_translate("Form", "Reference specific"))
		self.r_numReadsLineEdit_3.setText(_translate("Form", "--"))
		self.label_17.setText(_translate("Form", "Read number"))
		self.r_avQualLineEdit.setText(_translate("Form", "--"))
		self.label_18.setText(_translate("Form", "Average  quality"))
		self.r_estCoverageLineEdit.setText(_translate("Form", "--"))
		self.label_19.setText(_translate("Form", "Coverage"))
		self.readsFileButton.setText(_translate("Form", "Select"))
		self.qualityStatsButton.setText(_translate("Form", "Plot quality score distribution"))
		self.filterButton.setText(_translate("Form", "Filter"))
		self.qualityLineEdit.setText(_translate("Form", "30"))
		self.label_9.setText(_translate("Form", "Min. quality"))
	

	



if __name__ == "__main__":
	import sys
	app = QtWidgets.QApplication(sys.argv)
	Form = QtWidgets.QWidget()
	ui = Ui_Form()
	ui.setupUi(Form,installationDirectory)
	Form.show()
	sys.exit(app.exec_())