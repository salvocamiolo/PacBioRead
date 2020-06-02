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
		Form.resize(900, 518)
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
		self.label_6 = QtWidgets.QLabel(Form)
		self.label_6.setGeometry(QtCore.QRect(30, 220, 31, 31))
		self.label_6.setText("")
		self.label_6.setPixmap(QtGui.QPixmap(installationDirectory+"src/Images/1024px-Green_tick.png"))
		self.label_6.setObjectName("label_6")
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
		self.label_10 = QtWidgets.QLabel(Form)
		self.label_10.setGeometry(QtCore.QRect(30, 300, 31, 31))
		self.label_10.setText("")
		self.label_10.setPixmap(QtGui.QPixmap(installationDirectory+"src/Images/1024px-Green_tick.png"))
		self.label_10.setObjectName("label_10")
		self.readFilteringCheckBox = QtWidgets.QCheckBox(Form)
		self.readFilteringCheckBox.setGeometry(QtCore.QRect(70, 310, 211, 20))
		self.readFilteringCheckBox.setObjectName("readFilteringCheckBox")
		self.coverageLineEdit = QtWidgets.QLineEdit(Form)
		self.coverageLineEdit.setGeometry(QtCore.QRect(290, 310, 101, 21))
		self.coverageLineEdit.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.coverageLineEdit.setObjectName("coverageLineEdit")
		self.label_11 = QtWidgets.QLabel(Form)
		self.label_11.setGeometry(QtCore.QRect(290, 290, 101, 16))
		self.label_11.setObjectName("label_11")
		self.readFilteringCheckBox_2 = QtWidgets.QCheckBox(Form)
		self.readFilteringCheckBox_2.setGeometry(QtCore.QRect(70, 390, 211, 20))
		self.readFilteringCheckBox_2.setObjectName("readFilteringCheckBox_2")
		self.label_12 = QtWidgets.QLabel(Form)
		self.label_12.setGeometry(QtCore.QRect(30, 380, 31, 31))
		self.label_12.setText("")
		self.label_12.setPixmap(QtGui.QPixmap(installationDirectory+"src/Images/1024px-Green_tick.png"))
		self.label_12.setObjectName("label_12")
		self.label_13 = QtWidgets.QLabel(Form)
		self.label_13.setGeometry(QtCore.QRect(290, 370, 101, 16))
		self.label_13.setObjectName("label_13")
		self.windowSizeLineEdit = QtWidgets.QLineEdit(Form)
		self.windowSizeLineEdit.setGeometry(QtCore.QRect(290, 390, 101, 21))
		self.windowSizeLineEdit.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.windowSizeLineEdit.setObjectName("windowSizeLineEdit")
		self.windowStepLineEdit = QtWidgets.QLineEdit(Form)
		self.windowStepLineEdit.setGeometry(QtCore.QRect(430, 390, 101, 21))
		self.windowStepLineEdit.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.windowStepLineEdit.setObjectName("windowStepLineEdit")
		self.label_14 = QtWidgets.QLabel(Form)
		self.label_14.setGeometry(QtCore.QRect(430, 370, 101, 16))
		self.label_14.setObjectName("label_14")
		self.readFilteringCheckBox_3 = QtWidgets.QCheckBox(Form)
		self.readFilteringCheckBox_3.setGeometry(QtCore.QRect(70, 470, 211, 20))
		self.readFilteringCheckBox_3.setObjectName("readFilteringCheckBox_3")
		self.label_15 = QtWidgets.QLabel(Form)
		self.label_15.setGeometry(QtCore.QRect(30, 460, 31, 31))
		self.label_15.setText("")
		self.label_15.setPixmap(QtGui.QPixmap(installationDirectory+"src/Images/1024px-Green_tick.png"))
		self.label_15.setObjectName("label_15")
		self.runButton = QtWidgets.QPushButton(Form)
		self.runButton.setGeometry(QtCore.QRect(770, 450, 113, 32))
		self.runButton.setObjectName("runButton")
		self.exitButton = QtWidgets.QPushButton(Form)
		self.exitButton.setGeometry(QtCore.QRect(770, 480, 113, 32))
		self.exitButton.setObjectName("exitButton")
		self.plainTextEdit = QtWidgets.QTextEdit(Form)
		self.plainTextEdit.setGeometry(QtCore.QRect(570, 30, 311, 411))
		self.plainTextEdit.setObjectName("plainTextEdit")
		self.label_16 = QtWidgets.QLabel(Form)
		self.label_16.setGeometry(QtCore.QRect(570, 10, 111, 16))
		self.label_16.setObjectName("label_16")

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
		if os.path.isdir(outputFolder) == False:
			os.system("mkdir "+outputFolder)

		
		# Perform HQ read fragmentaiton
		inputFile = self.readsFileLineEdit.text()
		threshold = self.qualityLineEdit.text()
		minLen = self.lengthLineEdit.text()

		numSeq = 0
		outfile = open(outputFolder+"/masked.fasta","w")
		for seq_record in SeqIO.parse(inputFile,"fastq"):
			numSeq+=1
			if numSeq%1000 == 0:
				self.plainTextEdit.append(str(numSeq)+" analyzed....")
				self.plainTextEdit.repaint()
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
		self.readFilteringCheckBox.setText(_translate("Form", "Filter arganism specific reads"))
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



	

	



if __name__ == "__main__":
	import sys
	app = QtWidgets.QApplication(sys.argv)
	Form = QtWidgets.QWidget()
	ui = Ui_Form()
	ui.setupUi(Form,installationDirectory)
	Form.show()
	sys.exit(app.exec_())