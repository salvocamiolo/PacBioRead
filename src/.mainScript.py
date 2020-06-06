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



class Ui_Form(object):
	def setupUi(self, Form,installationDirectory):
		Form.setObjectName("Form")
		Form.resize(900, 651)
		self.label = QtWidgets.QLabel(Form)
		self.label.setGeometry(QtCore.QRect(10, 10, 101, 16))
		self.label.setObjectName("label")
		self.projectFolderLineEdit = QtWidgets.QLineEdit(Form)
		self.projectFolderLineEdit.setGeometry(QtCore.QRect(10, 30, 361, 21))
		self.projectFolderLineEdit.setObjectName("projectFolderLineEdit")
		self.projectFolderButton = QtWidgets.QPushButton(Form)
		self.projectFolderButton.setGeometry(QtCore.QRect(390, 27, 113, 30))
		self.projectFolderButton.setObjectName("projectFolderButton")
		self.readsFileButton = QtWidgets.QPushButton(Form)
		self.readsFileButton.setGeometry(QtCore.QRect(390, 78, 113, 30))
		self.readsFileButton.setObjectName("readsFileButton")
		self.label_2 = QtWidgets.QLabel(Form)
		self.label_2.setGeometry(QtCore.QRect(10, 60, 101, 16))
		self.label_2.setObjectName("label_2")
		self.readsFileLineEdit = QtWidgets.QLineEdit(Form)
		self.readsFileLineEdit.setGeometry(QtCore.QRect(10, 80, 361, 21))
		self.readsFileLineEdit.setObjectName("readsFileLineEdit")
		self.referenceButton = QtWidgets.QPushButton(Form)
		self.referenceButton.setGeometry(QtCore.QRect(390, 127, 113, 30))
		self.referenceButton.setObjectName("referenceButton")
		self.label_3 = QtWidgets.QLabel(Form)
		self.label_3.setGeometry(QtCore.QRect(10, 110, 101, 16))
		self.label_3.setObjectName("label_3")
		self.referenceLineEdit = QtWidgets.QLineEdit(Form)
		self.referenceLineEdit.setGeometry(QtCore.QRect(10, 130, 361, 21))
		self.referenceLineEdit.setObjectName("referenceLineEdit")
		self.label_4 = QtWidgets.QLabel(Form)
		self.label_4.setGeometry(QtCore.QRect(90, 170, 111, 16))
		self.label_4.setObjectName("label_4")
		self.label_5 = QtWidgets.QLabel(Form)
		self.label_5.setGeometry(QtCore.QRect(20, 170, 51, 16))
		self.label_5.setObjectName("label_5")
		self.hqFragmentationLabel = QtWidgets.QLabel(Form)
		self.hqFragmentationLabel.setGeometry(QtCore.QRect(30, 220, 31, 31))
		self.hqFragmentationLabel.setText("")
		self.hqFragmentationLabel.setPixmap(QtGui.QPixmap(installationDirectory+"/src/Images/1024px-Green_tick.png"))
		self.hqFragmentationLabel.setObjectName("hqFragmentationLabel")
		self.hqFragmentationCheckBox = QtWidgets.QCheckBox(Form)
		self.hqFragmentationCheckBox.setGeometry(QtCore.QRect(70, 230, 151, 20))
		self.hqFragmentationCheckBox.setObjectName("hqFragmentationCheckBox")
		self.label_7 = QtWidgets.QLabel(Form)
		self.label_7.setGeometry(QtCore.QRect(380, 170, 111, 16))
		self.label_7.setObjectName("label_7")
		self.label_8 = QtWidgets.QLabel(Form)
		self.label_8.setGeometry(QtCore.QRect(290, 210, 101, 16))
		self.label_8.setObjectName("label_8")
		self.label_9 = QtWidgets.QLabel(Form)
		self.label_9.setGeometry(QtCore.QRect(430, 210, 101, 16))
		self.label_9.setObjectName("label_9")
		self.lengthLineEdit = QtWidgets.QLineEdit(Form)
		self.lengthLineEdit.setGeometry(QtCore.QRect(290, 230, 101, 21))
		self.lengthLineEdit.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.lengthLineEdit.setObjectName("lengthLineEdit")
		self.qualityLineEdit = QtWidgets.QLineEdit(Form)
		self.qualityLineEdit.setGeometry(QtCore.QRect(430, 230, 101, 21))
		self.qualityLineEdit.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.qualityLineEdit.setObjectName("qualityLineEdit")
		self.referenceReadMappingLabel = QtWidgets.QLabel(Form)
		self.referenceReadMappingLabel.setGeometry(QtCore.QRect(30, 300, 31, 31))
		self.referenceReadMappingLabel.setText("")
		self.referenceReadMappingLabel.setPixmap(QtGui.QPixmap(installationDirectory+"/src/Images/1024px-Green_tick.png"))
		self.referenceReadMappingLabel.setObjectName("referenceReadMappingLabel")
		self.readsMappingCheckBox = QtWidgets.QCheckBox(Form)
		self.readsMappingCheckBox.setGeometry(QtCore.QRect(70, 310, 211, 20))
		self.readsMappingCheckBox.setObjectName("readsMappingCheckBox")
		self.coverageLineEdit = QtWidgets.QLineEdit(Form)
		self.coverageLineEdit.setGeometry(QtCore.QRect(290, 390, 101, 21))
		self.coverageLineEdit.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.coverageLineEdit.setObjectName("coverageLineEdit")
		self.label_11 = QtWidgets.QLabel(Form)
		self.label_11.setGeometry(QtCore.QRect(290, 370, 101, 16))
		self.label_11.setObjectName("label_11")
		self.readFilteringCheckBox_2 = QtWidgets.QCheckBox(Form)
		self.readFilteringCheckBox_2.setGeometry(QtCore.QRect(70, 470, 211, 20))
		self.readFilteringCheckBox_2.setObjectName("readFilteringCheckBox_2")
		self.denovoAssemblyLabel = QtWidgets.QLabel(Form)
		self.denovoAssemblyLabel.setGeometry(QtCore.QRect(30, 460, 31, 31))
		self.denovoAssemblyLabel.setText("")
		self.denovoAssemblyLabel.setPixmap(QtGui.QPixmap(installationDirectory+"/src/Images/1024px-Green_tick.png"))
		self.denovoAssemblyLabel.setObjectName("denovoAssemblyLabel")
		self.label_13 = QtWidgets.QLabel(Form)
		self.label_13.setGeometry(QtCore.QRect(290, 450, 101, 16))
		self.label_13.setObjectName("label_13")
		self.windowSizeLineEdit = QtWidgets.QLineEdit(Form)
		self.windowSizeLineEdit.setGeometry(QtCore.QRect(290, 470, 101, 21))
		self.windowSizeLineEdit.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.windowSizeLineEdit.setObjectName("windowSizeLineEdit")
		self.windowStepLineEdit = QtWidgets.QLineEdit(Form)
		self.windowStepLineEdit.setGeometry(QtCore.QRect(430, 470, 101, 21))
		self.windowStepLineEdit.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.windowStepLineEdit.setObjectName("windowStepLineEdit")
		self.label_14 = QtWidgets.QLabel(Form)
		self.label_14.setGeometry(QtCore.QRect(430, 450, 101, 16))
		self.label_14.setObjectName("label_14")
		self.readFilteringCheckBox_3 = QtWidgets.QCheckBox(Form)
		self.readFilteringCheckBox_3.setGeometry(QtCore.QRect(70, 550, 211, 20))
		self.readFilteringCheckBox_3.setObjectName("readFilteringCheckBox_3")
		self.consensusCorrectionLaebl = QtWidgets.QLabel(Form)
		self.consensusCorrectionLaebl.setGeometry(QtCore.QRect(30, 540, 31, 31))
		self.consensusCorrectionLaebl.setText("")
		self.consensusCorrectionLaebl.setPixmap(QtGui.QPixmap(installationDirectory+"/src/Images/1024px-Green_tick.png"))
		self.consensusCorrectionLaebl.setObjectName("consensusCorrectionLaebl")
		self.runButton = QtWidgets.QPushButton(Form)
		self.runButton.setGeometry(QtCore.QRect(650, 610, 113, 32))
		self.runButton.setObjectName("runButton")
		self.exitButton = QtWidgets.QPushButton(Form)
		self.exitButton.setGeometry(QtCore.QRect(770, 610, 113, 32))
		self.exitButton.setObjectName("exitButton")
		self.label_16 = QtWidgets.QLabel(Form)
		self.label_16.setGeometry(QtCore.QRect(570, 10, 111, 16))
		self.label_16.setObjectName("label_16")
		self.readFilteringCheckBox = QtWidgets.QCheckBox(Form)
		self.readFilteringCheckBox.setGeometry(QtCore.QRect(70, 390, 211, 20))
		self.readFilteringCheckBox.setObjectName("readFilteringCheckBox")
		self.specificReadsLabel = QtWidgets.QLabel(Form)
		self.specificReadsLabel.setGeometry(QtCore.QRect(30, 380, 31, 31))
		self.specificReadsLabel.setText("")
		self.specificReadsLabel.setPixmap(QtGui.QPixmap(installationDirectory+"/src/Images/1024px-Green_tick.png"))
		self.specificReadsLabel.setObjectName("specificReadsLabel")
		self.logTextEdit = QtWidgets.QTextEdit(Form)
		self.logTextEdit.setGeometry(QtCore.QRect(560, 30, 321, 561))
		self.logTextEdit.setObjectName("logTextEdit")

		self.retranslateUi(Form)
		QtCore.QMetaObject.connectSlotsByName(Form)


		self.projectFolderButton.clicked.connect(self.selectProjectFolder)
		self.readsFileButton.clicked.connect(self.selectReadFile)
		self.referenceButton.clicked.connect(self.selectReferenceFile)
		self.runButton.clicked.connect(self.runAssembly)
		self.exitButton.clicked.connect(self.exitProgram)


	def selectProjectFolder(self):
		folderName = QFileDialog.getExistingDirectory(None, "Select project folder","./")
		self.projectFolderLineEdit.setText(folderName)
		self.hqFragmentationLabel.setPixmap(QtGui.QPixmap(installationDirectory+"/src/Images/1024px-Red_tick.png"))
		self.referenceReadMappingLabel.setPixmap(QtGui.QPixmap(installationDirectory+"/src/Images/1024px-Red_tick.png"))
		self.specificReadsLabel.setPixmap(QtGui.QPixmap(installationDirectory+"/src/Images/1024px-Red_tick.png"))
		self.denovoAssemblyLabel.setPixmap(QtGui.QPixmap(installationDirectory+"/src/Images/1024px-Red_tick.png"))
		self.consensusCorrectionLaebl.setPixmap(QtGui.QPixmap(installationDirectory+"/src/Images/1024px-Red_tick.png"))

		if os.path.isdir(folderName) == False:
			os.system("mkdir "+folderName)

		else:
			if os.path.isfile(folderName+"/hq.fasta") == True:
				self.hqFragmentationLabel.setPixmap(QtGui.QPixmap(installationDirectory+"/src/Images/1024px-Green_tick.png"))
			
			if os.path.isfile(folderName+"/lastzOutput.txt") == True:
				self.referenceReadMappingLabel.setPixmap(QtGui.QPixmap(installationDirectory+"/src/Images/1024px-Green_tick.png"))
			
			if os.path.isfile(folderName+"/hq_specific.fasta") == True:
				self.specificReadsLabel.setPixmap(QtGui.QPixmap(installationDirectory+"/src/Images/1024px-Green_tick.png"))

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

		#Create project folder if it does not exist
		outputFolder = self.projectFolderLineEdit.text()

		# Perform HQ read fragmentaiton
		if self.hqFragmentationCheckBox.isChecked() == True:
			self.hqFragmentationLabel.setPixmap(QtGui.QPixmap(installationDirectory+"/src/Images/1024px-Yellow_tick.png"))
			self.logTextEdit.append("High quality read fragmentation started\n")
			self.logTextEdit.repaint()

			inputFile = self.readsFileLineEdit.text()
			threshold = self.qualityLineEdit.text()
			minLen = self.lengthLineEdit.text()

			numSeq = 0
			outfile = open(outputFolder+"/masked.fasta","w")
			for seq_record in SeqIO.parse(inputFile,"fastq"):
				numSeq+=1
				if numSeq%1000 == 0:
					self.logTextEdit.append(str(numSeq)+" analyzed....")
					self.logTextEdit.repaint()
					print(str(numSeq)+" analyzed....\n")
				sequence = str(seq_record.seq)
				quality = seq_record.letter_annotations["phred_quality"]
				maskedSeq = ""
				maskedQual = ""
				for a in range(len(quality)):
					if quality[a]>int(threshold):
						maskedSeq+=sequence[a]
						maskedQual+="G"
					else:
						outfile.write(">MaskedSeq_"+str(numSeq)+"\n"+maskedSeq+"\n")
						numSeq+=1
						maskedSeq=""
			outfile.close()
			outfile = open(outputFolder+"/hq.fasta","w")
			for seq_record in SeqIO.parse(outputFolder+"/masked.fasta","fasta"):
				if len(str(seq_record.seq))>int(minLen):
					SeqIO.write(seq_record,outfile,"fasta")

			os.system("rm "+outputFolder+"/masked.fasta")
			self.logTextEdit.append("\nHigh quality read fragmentation finished\n")
			self.logTextEdit.repaint()
			self.hqFragmentationLabel.setPixmap(QtGui.QPixmap(installationDirectory+"/src/Images/1024px-Green_tick.png"))
		


		#Perform read mapping on reference to extract organisms specific reads
		if self.readsMappingCheckBox.isChecked() == True:
			self.readsMappingCheckBox.setPixmap(QtGui.QPixmap(installationDirectory+"/src/Images/1024px-Yellow_tick.png"))
			self.logTextEdit.append("Read mapping on reference started\n")
			self.logTextEdit.repaint()
			if os.path.isfile(outputFolder+"/hq.fasta")==True:
				inputFile = outputFolder+"/hq.fasta"
				reference = self.referenceLineEdit.text()

				numSeq = 0
				os.system("rm -f "+outputFolder+"/lastzOutput.txt")
				os.system("touch "+outputFolder+"/lastzOutput.txt")
				for seq_record in SeqIO.parse(inputFile,"fasta"):
					outfile = open(outputFolder+"/tempFasta.fasta","w")
					SeqIO.write(seq_record,outfile,"fasta")
					outfile.close()
					os.system(installationDirectory+"/src/conda/bin/lastz "+outputFolder+"/tempFasta.fasta "+reference+"  --format=general  >> "+outputFolder+"/lastzOutput.txt")
					numSeq+=1
					print(numSeq)
				os.system("rm "+outputFolder+"/temp.fasta")
				self.logTextEdit.append("\nRead mapping on reference finished\n")
				self.logTextEdit.repaint()
				self.readsMappingCheckBox.setPixmap(QtGui.QPixmap(installationDirectory+"/src/Images/1024px-Green_tick.png"))
			else:
				msg = QMessageBox()
				msg.setIcon(QMessageBox.Warning)
				msg.setText("Please run the HQ fragmentation first")
				msg.setWindowTitle("Warning")
				msg.setDetailedText("You should perform the HQ fragmentation first. This will create a file named hq.fasta in your project folder, which is not there at the moment\n ")
				msg.setStandardButtons(QMessageBox.Ok)
				msg.exec_()
				return

		#Perform read mapping on reference to extract organisms specific reads
		if self.readFilteringCheckBox.isChecked() == True:
			if os.path.isfile(outputFolder+"/lastzOutput.txt") == True:
				if os.path.isfile(outputFolder+"/hq.fasta") == True:
					self.specificReadsLabel.setPixmap(QtGui.QPixmap(installationDirectory+"/src/Images/1024px-Yellow_tick.png"))
					inputFile = outputFolder+"/hq.fasta"
					lastzOutput = outputFolder+"/lastzOutput.txt"
					coverage = float(self.coverageLineEdit.text())/100

					self.logTextEdit.append("\nOrganisms specific read filtering started\n")
					self.logTextEdit.repaint()


					sequences = {}
					self.logTextEdit.append("\nLoading HQ reads in memory\n")
					self.logTextEdit.repaint()

					for seq_record in SeqIO.parse(inputFile,"fasta"):
						if not str(seq_record.id) in sequences:
							sequences[str(seq_record.id)] = seq_record

					self.logTextEdit.append("\nApplying filter\n")
					self.logTextEdit.repaint()
					outfile = open(outputFolder+"/hq_specific.fasta","w")
					infile = open(lastzOutput)
					numSeq = 0
					while True:
						line = infile.readline().rstrip()
						if not line:
							break
						if not "#score" in line:
							fields = line.split("\t")
							readName = fields[1]
							strand = fields[2]
							length = float(fields[3])
							start = float(fields[4])
							end = float(fields[5])

							if abs(end-start)/length> coverage:
								numSeq +=1
								if numSeq%1000 == 0:
									self.logTextEdit.append(str(numSeq)+" reads passed the filter")
									self.logTextEdit.repaint()
								SeqIO.write(sequences[readName],outfile,"fasta")
								print("Read %s comes from Reference!" %readName)
					self.logTextEdit.append(str(numSeq)+" reads passed the filter")
					self.logTextEdit.repaint()
					self.specificReadsLabel.setPixmap(QtGui.QPixmap(installationDirectory+"/src/Images/1024px-Green_tick.png"))
					infile.close()
					outfile.close()
				else:
					msg = QMessageBox()
					msg.setIcon(QMessageBox.Warning)
					msg.setText("Please run the HQ fragmentation first")
					msg.setWindowTitle("Warning")
					msg.setDetailedText("You should perform the HQ fragmentation first. This will create a file named hq.fasta in your project folder, which is not there at the moment\n ")
					msg.setStandardButtons(QMessageBox.Ok)
					msg.setStandardButtons(QMessageBox.Ok)
					msg.exec_()
					return

			else:
				msg = QMessageBox()
				msg.setIcon(QMessageBox.Warning)
				msg.setText("Please map the reads to the reference first.")
				msg.setWindowTitle("Warning")
				msg.setDetailedText("You should map the reads to the provided reference first. This will generate a file named lastzOutput.txt in your project folder\n ")
				msg.setStandardButtons(QMessageBox.Ok)
				msg.exec_()
				return
		

			
		#Perform reference guided de novo assembly
		if self.readFilteringCheckBox_2.isChecked() == True:
			reads = ""
			if os.path.isfile(outputFolder+"/hq_specific.fasta") == True:
				reads = outputFolder+"/hq_specific.fasta"
			if os.path.isfile(outputFolder+"/hq.fasta") == True:
				reads = outputFolder+"/hq.fasta"
			
			if reads == "":
				msg = QMessageBox()
				msg.setIcon(QMessageBox.Warning)
				msg.setText("Please run the HQ fragmentation first")
				msg.setWindowTitle("Warning")
				msg.setDetailedText("You should perform the HQ fragmentation first. This will create a file named hq.fasta in your project folder, which is not there at the moment\n ")
				msg.setStandardButtons(QMessageBox.Ok)
				msg.setStandardButtons(QMessageBox.Ok)
				msg.exec_()
				return

			else:
				refFile = self.referenceLineEdit.text()
				readsSeq = {}
				windowSize = int(self.windowSizeLineEdit.text())
				windowStep = int(self.windowStepLineEdit.text())

				self.logTextEdit.append("Reference guided de novo assembly started")
				self.logTextEdit.append("Loading reads in memory")

				self.logTextEdit.repaint()
				for seq_record in SeqIO.parse(reads,"fasta"):
					if not str(seq_record.id) in readsSeq:
						readsSeq[str(seq_record.id)] = str(seq_record.seq)

				reference = outputFolder+"/partReference.fasta"



				for seq_record in SeqIO.parse(refFile,"fasta"):
					refSeq = str(seq_record.seq)

				stage_a = open(outputFolder+"/preliminaryContigs.fasta","w")
				for a in range(0,len(refSeq),+windowStep):
					endPos = a+windowSize
					if endPos>len(refSeq):
						endPos=len(refSeq)

					self.logTextEdit.append("Assembling region "+str(a)+"-"+str(endPos))
					self.logTextEdit.repaint()
					#print("Assembling range %d-%d" %(a,a+windowSize))
					partSeq = refSeq[a:endPos]
					tempFasta = open(outputFolder+"/partReference.fasta","w")
					tempFasta.write(">partReference\n"+partSeq+"\n")
					tempFasta.close()

					self.logTextEdit.append("Aligning reads.... ")
					self.logTextEdit.repaint()

					with open(reads, "r") as fasta, open(reads+".fastq", "w") as fastq:
						for record in SeqIO.parse(fasta, "fasta"):
							record.letter_annotations["phred_quality"] = [40] * len(record)
							SeqIO.write(record, fastq, "fastq")
					os.system(installationDirectory+"/src/conda/bin/bowtie2-build "+reference+" bowtie2Ref")
					os.system(installationDirectory+"/src/conda/bin/bowtie2 -U "+reads+".fastq "+" "+" -x "+outputFolder+"/bowtie2Ref -S "+outputFolder+"/bowtie2Alignment.sam -p 8") #To add num threads
					os.system(installationDirectory+"/src/conda/bin/samtools view -F 4 -bS -h "+outputFolder+"/bowtie2Alignment.sam > "+outputFolder+"/bowtie2Mapped.bam")
					os.system(installationDirectory+"/src/conda/bin/samtools sort -o "+outputFolder+"/bowtie2Mapped_sorted.bam "+outputFolder+"/bowtie2Mapped.bam") 
					os.system(installationDirectory+"/src/conda/bin/samtools index "+outputFolder+"/bowtie2Mapped_sorted.bam")
					
					self.logTextEdit.append("Scanning alignment.... ")
					self.logTextEdit.repaint()
					readsToAssemble = set()
					for b in range(0,19500,+150):
						print("Analyzing range %d-%d" %(b,b+150))
						os.system(installationDirectory+"/src/conda/bin/samtools view "+outputFolder+"/bowtie2Mapped_sorted.bam partReference:"+str(b)+"-"+str(b+150)+" > "+outputFolder+"/localAlignment.sam")
						infile = open(outputFolder+"/localAlignment.sam")
						longestRead = ""
						longestReadLength = 0
						while True:
							line = infile.readline().rstrip()
							if not line:
								break
							fields = line.split("\t")
							if len(fields[7])>longestReadLength:
								longestReadLength = len(fields[7])
								longestRead = fields[0]
						readsToAssemble.add(longestRead)
						print(longestRead,longestReadLength)
						print("Selected %d sequences" %len(readsToAssemble))
						infile.close()

					outfile = open(outputFolder+"/toAssemble.fasta","w")
					numReadsToAssemble = 0
					for item in readsToAssemble:
						if not item == '':
							outfile.write(">"+item+"\n"+readsSeq[item]+"\n")
							numReadsToAssemble+=1
					outfile.close()
					print("Assembling %d reads with cap3" %numReadsToAssemble)
					self.logTextEdit.append("Assembling "+str(numReadsToAssemble)+" reads with cap3")
					self.logTextEdit.repaint()
					os.system(installationDirectory+"/src/conda/bin/cap3 "+outputFolder+"/toAssemble.fasta >null 2>&1")
					maxScaffoldLength = 0
					longestContig = ""
					for seq_record in SeqIO.parse(outputFolder+"/toAssemble.fasta.cap.contigs","fasta"):
						if len(str(seq_record.seq)) > maxScaffoldLength:
							maxScaffoldLength = len(str(seq_record.seq))
							longestContig = str(seq_record.seq)

					self.logTextEdit.append("Scaffold size: "+str(maxScaffoldLength))
					self.logTextEdit.repaint()
					stage_a.write(">Range_"+str(a)+"_"+str(endPos)+"\n"+longestContig+"\n")

				stage_a.close()
				os.system("rm "+outputFolder+"/partReference.fasta* "+outputFolder+"/outputBlast.txt "+outputFolder+"/toAssemble*")

				#Assembling all the contigs with cap3 to obtain extended contigs
				self.logTextEdit.append("\nJoining contigs.... ")
				self.logTextEdit.repaint()
				os.system(installationDirectory+"/src/conda/bin/cap3 "+outputFolder+"/preliminaryContigs.fasta")

				
				
				#Scaffolding all the extended contigs in the previous step with ragout
				self.logTextEdit.append("\nScaffolding.... ")
				self.logTextEdit.repaint()
				ragoutRecepie = open(outputFolder+"/ragout_recepie.rcp","w")
				ragoutRecepie.write(".references = reference\n.target = scaffolds\n\nreference.fasta = "+refFile+"\nscaffolds.fasta = "+outputFolder+"/preliminaryContigs.fasta.cap.contigs")
				ragoutRecepie.close()
				os.system(installationDirectory+"/src/conda2/bin/ragout -o "+outputFolder+"/ragoutOutput "+outputFolder+"/ragout_recepie.rcp")
				os.system("cp "+outputFolder+"/preliminaryContigs.fasta "+outputFolder+"/stage_a.fasta")
				os.system("cp "+outputFolder+"/ragoutOutput/scaffolds_scaffolds.fasta "+outputFolder+"/stage_b.fasta")
				os.system("rm -rf bowtie2* null local*")

				#Check the present of N and if present close the gaps with lr_gapcloser
				for seq_record in SeqIO.parse(outputFolder+"/stage_b.fasta","fasta"):
					scaffoldSseq = str(seq_record.seq)
					position = -1
					while position < len(str(seq_record.seq))-1:
						position+=1
						if scaffoldSseq[position] == "N" or scaffoldSseq[position] == "n":
							gapStart = position
							while scaffoldSseq[position] == "N" or scaffoldSseq[position] == 'n':
								position+=1
							gapEnd = position
							print("Found gap between position %d and %d " %(gapStart, gapEnd))
							bitToJoin = open(outputFolder+"/firstBit.fasta","w")
							bitToJoin.write(">firstBit\n"+scaffoldSseq[gapStart-2000:gapStart]+"\n")
							bitToJoin.close()
							bitToJoin = open(outputFolder+"/secondBit.fasta","w")
							bitToJoin.write(">secondBit\n"+scaffoldSseq[gapEnd:gapEnd+500]+"\n")
							bitToJoin.close()
							#converting original fastq file into fasta file
							self.logTextEdit.append("Closing gap.... ")
							self.logTextEdit.repaint()
							os.system(installationDirectory+"/src/conda/bin/python "+installationDirectory+"/src/scripts/lr_gapCloser.py -p " \
								+installationDirectory+" -i "+self.readsFileLineEdit.text()+" -s "+outputFolder+"/firstBit.fasta - e "+ \
									outputFolder+"/secondBit.fasta -x "+outputFolder+" -o "+"gap_"+str(gapStart)+"_"+str(gapEnd))
							os.system("cat "+outputFolder+"/gap_"+str(gapStart)+"_"+str(gapEnd)+" >> " \
								+outputFolder+"/preliminaryContigs.fasta.cap.contigs")

							os.system(installationDirectory+"/src/conda2/bin/ragout -o "+outputFolder+"/ragoutOutput_"+str(gapStart)+"_"+str(gapEnd)+" "+outputFolder+"/ragout_recepie.rcp")
							
							


					

	





			

			

			




			

		





			



	def retranslateUi(self, Form):
		_translate = QtCore.QCoreApplication.translate
		Form.setWindowTitle(_translate("Form", "Form"))
		self.label.setText(_translate("Form", "Project folder"))
		self.projectFolderButton.setText(_translate("Form", "Select"))
		self.readsFileButton.setText(_translate("Form", "Select"))
		self.label_2.setText(_translate("Form", "Reads file"))
		self.referenceButton.setText(_translate("Form", "Select"))
		self.label_3.setText(_translate("Form", "Reference file"))
		self.label_4.setText(_translate("Form", "Tasks to perform"))
		self.label_5.setText(_translate("Form", "Status"))
		self.hqFragmentationCheckBox.setText(_translate("Form", "HQ fragmentation"))
		self.label_7.setText(_translate("Form", "Parameters"))
		self.label_8.setText(_translate("Form", "Min. read length"))
		self.label_9.setText(_translate("Form", "Min. quality"))
		self.lengthLineEdit.setText(_translate("Form", "150"))
		self.qualityLineEdit.setText(_translate("Form", "30"))
		self.readsMappingCheckBox.setText(_translate("Form", "Map reads to reference"))
		self.coverageLineEdit.setText(_translate("Form", "90"))
		self.label_11.setText(_translate("Form", "Coverage (%)"))
		self.readFilteringCheckBox_2.setText(_translate("Form", "De novo assembly"))
		self.label_13.setText(_translate("Form", "Window size"))
		self.windowSizeLineEdit.setText(_translate("Form", "20000"))
		self.windowStepLineEdit.setText(_translate("Form", "10000"))
		self.label_14.setText(_translate("Form", "Window step"))
		self.readFilteringCheckBox_3.setText(_translate("Form", "Consensus correction"))
		self.runButton.setText(_translate("Form", "Run"))
		self.exitButton.setText(_translate("Form", "Exit"))
		self.label_16.setText(_translate("Form", "Log "))
		self.readFilteringCheckBox.setText(_translate("Form", "Select reference specific reads"))



	

	



if __name__ == "__main__":
	import sys
	app = QtWidgets.QApplication(sys.argv)
	Form = QtWidgets.QWidget()
	ui = Ui_Form()
	ui.setupUi(Form,installationDirectory)
	Form.show()
	sys.exit(app.exec_())