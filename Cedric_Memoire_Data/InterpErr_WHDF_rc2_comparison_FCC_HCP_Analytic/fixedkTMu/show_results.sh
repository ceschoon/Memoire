#! /bin/bash

outFile="results.dat"
echo "" > $outFile

echo "" >> $outFile
echo "Thermodynamic parameters" >> $outFile
cat dataFCC/kTMuSolid_kT\=* | grep kT >> $outFile
cat dataFCC/kTMuSolid_kT\=* | grep mu >> $outFile

echo "" >> $outFile
echo "FCC results" >> $outFile
cat dataFCC/kTMuSolid_kT\=* | grep freeEnergy >> $outFile
cat dataFCC/kTMuSolid_kT\=* | grep density >> $outFile

echo "" >> $outFile
echo "HCP results" >> $outFile
cat dataHCP/kTMuSolid_kT\=* | grep freeEnergy >> $outFile
cat dataHCP/kTMuSolid_kT\=* | grep density >> $outFile

cat $outFile
