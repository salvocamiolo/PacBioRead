# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file './src/GUI/GUI2.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Form(object):
    def setupUi(self, Form):
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

