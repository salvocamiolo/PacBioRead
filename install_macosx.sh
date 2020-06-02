installationPath=$(pwd)
echo "#!"$installationPath"/src/conda/bin/python" >PacBioReads
echo "installationDirectory = \""$installationPath"/\"" >> PacBioReads
echo " " >> PacBioReads
cat PacBioReads ./src/.mainScript >temp ; mv temp PacBioReads
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