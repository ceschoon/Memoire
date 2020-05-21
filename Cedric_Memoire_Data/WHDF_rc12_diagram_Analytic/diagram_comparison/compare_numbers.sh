#! /bin/bash

outFile="results.dat"
echo "" > $outFile

echo "" >> $outFile
echo "Numeric Weights:" >> $outFile
cat ../coexCurves/coexCurveFS.dat | head -n 5 >> $outFile

echo "" >> $outFile
echo "Analytic Weights:" >> $outFile
cat ../coexCurves/coexCurveFS.dat | head -n 5 >> $outFile

cat $outFile
