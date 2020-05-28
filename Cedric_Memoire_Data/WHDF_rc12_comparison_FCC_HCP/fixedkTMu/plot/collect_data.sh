#! /bin/bash

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

paste datakT.dat dataMu.dat dataFCC.dat dataHCP.dat > data.dat

rm datakT.dat
rm dataMu.dat
rm dataFCC.dat
rm dataHCP.dat

