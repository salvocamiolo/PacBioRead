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
		Form.resize(846, 520)
		self.label = QtWidgets.QLabel(Form)
		self.label.setGeometry(QtCore.QRect(10, 10, 101, 16))
		self.label.setObjectName("label")
		self.projectFolderLineEdit = QtWidgets.QLineEdit(Form)
		self.projectFolderLineEdit.setGeometry(QtCore.QRect(10, 30, 361, 21))
		self.projectFolderLineEdit.setObjectName("projectFolderLineEdit")
		self.projectFolderButton = QtWidgets.QPushButton(Form)
		self.projectFolderButton.setGeometry(QtCore.QRect(390, 30, 113, 21))
		self.projectFolderButton.setObjectName("projectFolderButton")
		self.readsFileButton = QtWidgets.QPushButton(Form)
		self.readsFileButton.setGeometry(QtCore.QRect(390, 80, 113, 21))
		self.readsFileButton.setObjectName("readsFileButton")
		self.label_2 = QtWidgets.QLabel(Form)
		self.label_2.setGeometry(QtCore.QRect(10, 60, 101, 16))
		self.label_2.setObjectName("label_2")
		self.readsFileLineEdit = QtWidgets.QLineEdit(Form)
		self.readsFileLineEdit.setGeometry(QtCore.QRect(10, 80, 361, 21))
		self.readsFileLineEdit.setObjectName("readsFileLineEdit")
		self.referenceButton = QtWidgets.QPushButton(Form)
		self.referenceButton.setGeometry(QtCore.QRect(390, 130, 113, 21))
		self.referenceButton.setObjectName("referenceButton")
		self.label_3 = QtWidgets.QLabel(Form)
		self.label_3.setGeometry(QtCore.QRect(10, 110, 101, 16))
		self.label_3.setObjectName("label_3")
		self.referenceLineEdit = QtWidgets.QLineEdit(Form)
		self.referenceLineEdit.setGeometry(QtCore.QRect(10, 130, 361, 21))
		self.referenceLineEdit.setObjectName("referenceLineEdit")
		self.label_9 = QtWidgets.QLabel(Form)
		self.label_9.setGeometry(QtCore.QRect(10, 230, 101, 16))
		self.label_9.setObjectName("label_9")
		self.qualityLineEdit = QtWidgets.QLineEdit(Form)
		self.qualityLineEdit.setGeometry(QtCore.QRect(10, 250, 101, 21))
		self.qualityLineEdit.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.qualityLineEdit.setObjectName("qualityLineEdit")
		self.label_13 = QtWidgets.QLabel(Form)
		self.label_13.setGeometry(QtCore.QRect(10, 170, 101, 16))
		self.label_13.setObjectName("label_13")
		self.windowSizeLineEdit = QtWidgets.QLineEdit(Form)
		self.windowSizeLineEdit.setGeometry(QtCore.QRect(10, 190, 101, 21))
		self.windowSizeLineEdit.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.windowSizeLineEdit.setObjectName("windowSizeLineEdit")
		self.windowStepLineEdit = QtWidgets.QLineEdit(Form)
		self.windowStepLineEdit.setGeometry(QtCore.QRect(140, 190, 101, 21))
		self.windowStepLineEdit.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.windowStepLineEdit.setObjectName("windowStepLineEdit")
		self.label_14 = QtWidgets.QLabel(Form)
		self.label_14.setGeometry(QtCore.QRect(140, 170, 101, 16))
		self.label_14.setObjectName("label_14")
		self.runButton = QtWidgets.QPushButton(Form)
		self.runButton.setGeometry(QtCore.QRect(600, 490, 113, 21))
		self.runButton.setObjectName("runButton")
		self.exitButton = QtWidgets.QPushButton(Form)
		self.exitButton.setGeometry(QtCore.QRect(720, 490, 113, 21))
		self.exitButton.setObjectName("exitButton")
		self.label_16 = QtWidgets.QLabel(Form)
		self.label_16.setGeometry(QtCore.QRect(20, 290, 111, 16))
		self.label_16.setObjectName("label_16")
		self.logTextEdit = QtWidgets.QTextEdit(Form)
		self.logTextEdit.setGeometry(QtCore.QRect(10, 310, 821, 161))
		self.logTextEdit.setObjectName("logTextEdit")
		self.qualityStatsButton = QtWidgets.QPushButton(Form)
		self.qualityStatsButton.setGeometry(QtCore.QRect(140, 250, 101, 21))
		self.qualityStatsButton.setObjectName("qualityStatsButton")
		self.label_17 = QtWidgets.QLabel(Form)
		self.label_17.setGeometry(QtCore.QRect(540, 10, 111, 16))
		self.label_17.setObjectName("label_17")
		self.textEdit = QtWidgets.QTextEdit(Form)
		self.textEdit.setGeometry(QtCore.QRect(540, 30, 291, 241))
		font = QtGui.QFont()
		font.setPointSize(8)
		self.textEdit.setFont(font)
		self.textEdit.setObjectName("textEdit")
		self.numThreadsLineEdit = QtWidgets.QLineEdit(Form)
		self.numThreadsLineEdit.setGeometry(QtCore.QRect(270, 190, 101, 21))
		self.numThreadsLineEdit.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.numThreadsLineEdit.setObjectName("numThreadsLineEdit")
		self.numThreadsLabel = QtWidgets.QLabel(Form)
		self.numThreadsLabel.setGeometry(QtCore.QRect(270, 170, 101, 16))
		self.numThreadsLabel.setObjectName("numThreadsLabel")

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

		#Convert the input fastq file in fast format
		outputFolder = self.projectFolderLineEdit.text()
		os.system(installationDirectory+"/src/conda/bin/fq2fa "+self.readsFileLineEdit.text()+" "+outputFolder+'/originalReads.fasta')


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
			self.referenceReadMappingLabel.setPixmap(QtGui.QPixmap(installationDirectory+"/src/Images/1024px-Yellow_tick.png"))
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
				self.referenceReadMappingLabel.setPixmap(QtGui.QPixmap(installationDirectory+"/src/Images/1024px-Green_tick.png"))
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
			self.denovoAssemblyLabel.setPixmap(QtGui.QPixmap(installationDirectory+"/src/Images/1024px-Yellow_tick.png"))
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
				#for seq_record in SeqIO.parse(self.readsFileLineEdit.text(),"fastq"):
					if not str(seq_record.id) in readsSeq:
						readsSeq[str(seq_record.id)] = str(seq_record.seq)

				reference = outputFolder+"/partReference.fasta"



				for seq_record in SeqIO.parse(refFile,"fasta"):
					refSeq = str(seq_record.seq)

				stage_a = open(outputFolder+"/stage_a.fasta","w")

				with open(reads, "r") as fasta, open(reads+".fastq", "w") as fastq:
					for record in SeqIO.parse(fasta, "fasta"):
						record.letter_annotations["phred_quality"] = [40] * len(record)
						SeqIO.write(record, fastq, "fastq")
				
				for a in range(0,len(refSeq),+windowStep):
				#for a in range(1):
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

					


					os.system(installationDirectory+"/src/conda/bin/minimap2 "+outputFolder+"/partReference.fasta "+reads+".fastq > "+outputFolder+"/outputMinimap")

					os.system("awk '(($4-$3)/$2)>0.80' "+outputFolder+"/outputMinimap | sort -k2rn,2rn >  "+outputFolder+"/outputMinimap_filtered ")

					self.logTextEdit.append("Scanning alignment.... ")
					self.logTextEdit.repaint()
					readsToAssemble = set()
					numAttempt = 0
					maxScaffoldLength = 0

					
					while float(maxScaffoldLength) < float(windowSize)*0.9:
						numAttempt +=1
						if numAttempt == 5:
							break
						
						for b in range(0,19500,+150):
							tfile = open(outputFolder+"/outputMinimap_filtered")
							print("Analyzing range %d-%d" %(b,b+150))
							
							refPos = -1
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

						print("Assembling %d reads with cap3" %numReadsToAssemble)
						self.logTextEdit.append("Assembling "+str(numReadsToAssemble)+" reads")
						self.logTextEdit.repaint()
						#os.system(installationDirectory+"/src/conda/bin/cap3 "+outputFolder+"/toAssemble.fasta >null 2>&1")
						os.system(installationDirectory+"/src/conda/bin/art_illumina -i "+outputFolder+"/toAssemble.fasta -l 150 -f 30 -ss HS25 -o "+outputFolder+"/simulatedReads -p -m 500 -s 50")
						toAssembleFile = open(outputFolder+"/allSimulated.fasta","w")
						os.system(installationDirectory+"/src/conda/bin/fq2fa --merge "+outputFolder+"/simulatedReads1.fq "+outputFolder+"/simulatedReads2.fq "+outputFolder+"/allSimulated.fasta")
						os.system("rm -rf "+outputFolder+"/outputIdba/")
						os.system(installationDirectory+"/src/conda/bin/idba_hybrid --reference "+outputFolder+"/partReference.fasta -r "+outputFolder+"/allSimulated.fasta --num_threads 8 -o "+outputFolder+"/outputIdba")
						maxScaffoldLength = 0
						longestContig = ""
						
						for seq_record in SeqIO.parse(outputFolder+"//outputIdba/scaffold.fa","fasta"):
							if len(str(seq_record.seq)) > maxScaffoldLength:
								maxScaffoldLength = len(str(seq_record.seq))
								longestContig = str(seq_record.seq)

						self.logTextEdit.append("Scaffold size: "+str(maxScaffoldLength))
						self.logTextEdit.repaint()
					stage_a.write(">Range_"+str(a)+"_"+str(endPos)+"\n"+longestContig+"\n")

				stage_a.close()

				#Joining contigs
				self.logTextEdit.append("\nJoining contigs.... ")
				self.logTextEdit.repaint()
				os.system(installationDirectory+"/src/conda/bin/python "+installationDirectory+"/src/scripts/contigsJoiner.py -c "+outputFolder+"/stage_a.fasta -r "+refFile+" -p "+installationDirectory+" -o "+outputFolder)

				#Check final number of scaffold and attempt gap closure if > 1
				finalScaffols = SeqIO.to_dict(SeqIO.parse(outputFolder+"/stage_b.fasta"))

				if len(finalScaffols)>1:
					self.logTextEdit.append("\nAttempting gap closure.... ")
					self.logTextEdit.repaint()
					os.system(installationDirectory+"/src/conda/bin/python "+installationDirectory+"/src/scripts/oneReadContigsJoiner.py \
						-p "+installationDirectory+ " -c "+outputFolder+"/stage_b.fasta -r "+refFile+" -x " + \
							outputFolder+" -s "+ self.readsFileLineEdit.text()+" -o stage_c.fasta")

		#Align original reads and create consensus
		#if self.readFilteringCheckBox_3.isChecked() == True:




	def retranslateUi(self, Form):
		_translate = QtCore.QCoreApplication.translate
		Form.setWindowTitle(_translate("Form", "Form"))
		self.label.setText(_translate("Form", "Project folder"))
		self.projectFolderButton.setText(_translate("Form", "Select"))
		self.readsFileButton.setText(_translate("Form", "Select"))
		self.label_2.setText(_translate("Form", "Reads file"))
		self.referenceButton.setText(_translate("Form", "Select"))
		self.label_3.setText(_translate("Form", "Reference file"))
		self.label_9.setText(_translate("Form", "Min. quality"))
		self.qualityLineEdit.setText(_translate("Form", "30"))
		self.label_13.setText(_translate("Form", "Window size"))
		self.windowSizeLineEdit.setText(_translate("Form", "20000"))
		self.windowStepLineEdit.setText(_translate("Form", "10000"))
		self.label_14.setText(_translate("Form", "Window step"))
		self.runButton.setText(_translate("Form", "Run"))
		self.exitButton.setText(_translate("Form", "Exit"))
		self.label_16.setText(_translate("Form", "Log "))
		self.qualityStatsButton.setText(_translate("Form", "Quality stats"))
		self.label_17.setText(_translate("Form", "Quality stats"))
		self.textEdit.setHtml(_translate("Form", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Sans Serif\'; font-size:8pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'.AppleSystemUIFont\'; font-size:10pt;\">Original reads</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'.AppleSystemUIFont\'; font-size:10pt;\">Number:   --</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'.AppleSystemUIFont\'; font-size:10pt;\">Average quality:  --</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'.AppleSystemUIFont\'; font-size:10pt;\">Bases (Q&gt;30):  --</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'.AppleSystemUIFont\'; font-size:10pt;\">Bases (Q&gt;20):  --</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'.AppleSystemUIFont\'; font-size:10pt;\">Bases (Q&gt;20):  --</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'.AppleSystemUIFont\'; font-size:10pt;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'.AppleSystemUIFont\'; font-size:10pt;\">Reference homologous</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'.AppleSystemUIFont\'; font-size:10pt;\">Number:  --</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'.AppleSystemUIFont\'; font-size:10pt;\">Average quality:  --</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'.AppleSystemUIFont\'; font-size:10pt;\">Bases (Q&gt;30):  --</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'.AppleSystemUIFont\'; font-size:10pt;\">Bases (Q&gt;20):  --</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'.AppleSystemUIFont\'; font-size:10pt;\">Bases (Q&gt;20):  --</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'.AppleSystemUIFont\'; font-size:10pt;\"><br /></p></body></html>"))

		self.numThreadsLineEdit.setText(_translate("Form", "8"))
		self.numThreadsLabel.setText(_translate("Form", "Num. threads"))
	

	



if __name__ == "__main__":
	import sys
	app = QtWidgets.QApplication(sys.argv)
	Form = QtWidgets.QWidget()
	ui = Ui_Form()
	ui.setupUi(Form,installationDirectory)
	Form.show()
	sys.exit(app.exec_())