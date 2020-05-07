#! /bin/bash

mkdir -p coexCurves
mkdir -p coexCurves/logs
mkdir -p critical_point
mkdir -p dataFS

### merge coexistence data

cp -r dataFS1/* dataFS/
cp -r dataFS2/* dataFS/

### generate coexistence curves

cd coexCurves

../../../Cedric_Memoire_Code/coexCurveVL/coexCurveVL \
	../input_files/input_coexCurveVL.dat
mv log.dat logs/log_coexCurveVL.dat

../../../Cedric_Memoire_Code/coexCurveFS/coexCurveFS \
	../input_files/input_coexCurveFS.dat
mv log.dat logs/log_coexCurveFS.dat

cd ..

### draw phase diagram

cd diagram
./plot_diagram.sh
cd ..

### compute critical point

cd critical_point
../../../Cedric_Memoire_Code/criticalPoint/criticalPoint \
	../input_files/input_criticalPoint.dat
mv log_kTc.dat log_criticalPoint.dat
cd ..
