#! /bin/bash

mkdir -p coexCurves_FCC
mkdir -p coexCurves_HCP
mkdir -p coexCurves_HCP/logs
mkdir -p dataVS
mkdir -p dataLS

### merge coexistence data (HCP)

cp -r dataVS1/* dataVS/

cp -r dataLS1/* dataLS/
cp -r dataLS2/* dataLS/

### generate coexistence curves (HCP)

cd coexCurves_HCP

../../../Cedric_Memoire_Code/coexCurveFS/coexCurveFS \
	../input_files/input_coexCurveVS.dat
mv log.dat logs/log_coexCurveVS.dat

../../../Cedric_Memoire_Code/coexCurveFS/coexCurveFS \
	../input_files/input_coexCurveLS.dat
mv log.dat logs/log_coexCurveLS.dat

cd ..

### get FCC coexCurves

cp -r ../WHDF_rc2_diagram/coexCurves/* ./coexCurves_FCC/

### plot diagram

cd diagram
./plot_diagram.sh
cd ..