#! /bin/bash

# save original data 

zip -r ../dataFCC.zip ../dataFCC
zip -r ../dataHCP.zip ../dataHCP



# remove unsuccessful computations

for file in $(ls ../dataFCC/ | grep kTMu)
do
	success_str=$(cat ../dataFCC/$file | grep success)
	unsuccessful_str="success = 0"
	
	if [ "$success_str" = "$unsuccessful_str" ]
	then
		rm -v ../dataFCC/$file
		rm -v ../dataHCP/$file
	fi
done

for file in $(ls ../dataHCP/ | grep kTMu)
do
	success_str=$(cat ../dataHCP/$file | grep success)
	unsuccessful_str="success = 0"
	
	if [ "$success_str" = "$unsuccessful_str" ]
	then
		rm -v ../dataFCC/$file
		rm -v ../dataHCP/$file
	fi
done



# collect data in seperate files for each property 

outFile="datakT.dat"
echo "# temperatures for WHDF rc=1.2" > $outFile
echo "" >> $outFile
cat ../dataFCC/kTMuSolid_kT\=* | grep kT >> $outFile


outFile="dataMu.dat"
echo "# chemical potentials for WHDF rc=1.2" > $outFile
echo "" >> $outFile
cat ../dataFCC/kTMuSolid_kT\=* | grep mu >> $outFile


outFile="dataFCC.dat"
echo "# FCC free energies for WHDF rc=1.2" > $outFile
echo "" >> $outFile
cat ../dataFCC/kTMuSolid_kT\=* | grep freeEnergy >> $outFile


outFile="dataHCP.dat"
echo "# HCP free energies for WHDF rc=1.2" > $outFile
echo "" >> $outFile
cat ../dataHCP/kTMuSolid_kT\=* | grep freeEnergy >> $outFile



# group in one file the data of each property

paste datakT.dat dataMu.dat dataFCC.dat dataHCP.dat > data.dat

rm datakT.dat
rm dataMu.dat
rm dataFCC.dat
rm dataHCP.dat

