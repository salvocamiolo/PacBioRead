installationPath=$(pwd)
echo "#!"$installationPath"/src/conda/bin/python" >PacBioReads
echo "installationDirectory = \""$installationPath"/\"" >> PacBioReads
echo " " >> PacBioReads
cat PacBioReads ./src/.mainScript.py >temp ; mv temp PacBioReads
chmod +x PacBioReads

installationPath=$(pwd)
echo "#!"$installationPath"/src/conda/bin/python" >PacBioReads_commmandLine
echo "installationDirectory = \""$installationPath"/\"" >> PacBioReads_commmandLine
echo " " >> PacBioReads_commmandLine
cat PacBioReads_commmandLine ./src/.clScript.py >temp ; mv temp PacBioReads_commmandLine
chmod +x PacBioReads_commmandLine

if  test -f "./src/conda/bin/conda"; then
	echo "Miniconda3 already installed"
else


	cd src
	bash Miniconda3-latest-Linux-x86_64.sh -b -p ./conda
	cd ../
	if  test -f "./src/conda/bin/conda"; then
	echo "Miniconda3 successfully installed"
	./src/conda/bin/conda install -c bioconda -y  lastz=1.0.4
	./src/conda/bin/conda install -c bioconda -y  cap3
	./src/conda/bin/conda install -c bioconda -y  blast=2.9.0
	./src/conda/bin/conda install -c bioconda -y  samtools=1.3.1
	./src/conda/bin/conda install -c anaconda -y pyqt=5.9.2
	./src/conda/bin/conda install -c bioconda -y  biopython=1.76
	./src/conda/bin/conda install -c bioconda -y bowtie2
	./src/conda/bin/conda install -c bioconda -y art
	./src/conda/bin/conda install -c bioconda -y idba
	./src/conda/bin/conda install -c bioconda -y minimap2
	./src/conda/bin/conda install -c bioconda -y varscan
	./src/conda/bin/conda install -c bioconda -y tabix
	./src/conda/bin/conda install -c bioconda -y bcftools
	./src/conda/bin/conda install -c bioconda -y  jellyfish=2.2.10
	./src/conda/bin/conda install -c bioconda -y  raven
	
	fi
fi




