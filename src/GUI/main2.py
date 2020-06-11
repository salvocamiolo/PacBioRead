# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'GUI2.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(846, 536)
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
        self.runButton.setGeometry(QtCore.QRect(600, 490, 113, 32))
        self.runButton.setObjectName("runButton")
        self.exitButton = QtWidgets.QPushButton(Form)
        self.exitButton.setGeometry(QtCore.QRect(720, 490, 113, 32))
        self.exitButton.setObjectName("exitButton")
        self.label_16 = QtWidgets.QLabel(Form)
        self.label_16.setGeometry(QtCore.QRect(20, 290, 111, 16))
        self.label_16.setObjectName("label_16")
        self.logTextEdit = QtWidgets.QTextEdit(Form)
        self.logTextEdit.setGeometry(QtCore.QRect(10, 310, 821, 161))
        self.logTextEdit.setObjectName("logTextEdit")
        self.qualityStatsButton = QtWidgets.QPushButton(Form)
        self.qualityStatsButton.setGeometry(QtCore.QRect(136, 247, 111, 30))
        self.qualityStatsButton.setObjectName("qualityStatsButton")
        self.label_17 = QtWidgets.QLabel(Form)
        self.label_17.setGeometry(QtCore.QRect(540, 10, 111, 16))
        self.label_17.setObjectName("label_17")
        self.textEdit = QtWidgets.QTextEdit(Form)
        self.textEdit.setGeometry(QtCore.QRect(540, 30, 291, 241))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.textEdit.setFont(font)
        self.textEdit.setObjectName("textEdit")

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

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
"</style></head><body style=\" font-family:\'.AppleSystemUIFont\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:13pt;\">Original reads</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:13pt;\">Number:   --</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:13pt;\">Average quality:  --</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:13pt;\">Bases (Q&gt;30):  --</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:13pt;\">Bases (Q&gt;20):  --</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:13pt;\">Bases (Q&gt;20):  --</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:13pt;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:13pt;\">Reference homologue</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:13pt;\">Number:  --</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:13pt;\">Average quality:  --</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:13pt;\">Bases (Q&gt;30):  --</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:13pt;\">Bases (Q&gt;20):  --</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:13pt;\">Bases (Q&gt;20):  --</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:13pt;\"><br /></p></body></html>"))

