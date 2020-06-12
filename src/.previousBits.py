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