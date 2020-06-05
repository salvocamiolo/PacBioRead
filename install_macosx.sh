installationPath=$(pwd)
echo "#!"$installationPath"/src/conda/bin/python" >PacBioReads
echo "installationDirectory = \""$installationPath"/\"" >> PacBioReads
echo " " >> PacBioReads
cat PacBioReads ./src/.mainScript.py >temp ; mv temp PacBioReads
chmod +x PacBioReads

if  test -f "./src/conda/bin/conda"; then
	echo "Miniconda3 already installed"
else


	cd src
	bash Miniconda3-latest-MacOSX-x86_64.sh -b -p ./conda
	cd ../
	if  test -f "./src/conda/bin/conda"; then
	echo "Miniconda3 successfully installed"
	fi
fi


./src/conda/bin/conda install -c bioconda -y  lastz=1.0.4
./src/conda/bin/conda install -c bioconda -y  cap3
./src/conda/bin/conda install -c bioconda -y  blast=2.9.0
./src/conda/bin/conda install -c bioconda -y  samtools=1.3.1
./src/conda/bin/conda install -c anaconda -y pyqt=5.9.2
./src/conda/bin/conda install -c bioconda -y  biopython=1.76
./src/conda/bin/conda install -c bioconda bowtie2