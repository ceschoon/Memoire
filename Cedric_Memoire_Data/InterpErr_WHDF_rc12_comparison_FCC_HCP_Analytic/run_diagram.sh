#! /bin/bash

mkdir -p coexCurves_FCC
mkdir -p coexCurves_HCP
mkdir -p coexCurves_HCP/logs
mkdir -p dataFS

### merge coexistence data (HCP)

cp -r dataFS1/* dataFS/
cp -r dataFS2/* dataFS/
cp -r dataFS3/* dataFS/

### generate coexistence curves (HCP)

cd coexCurves_HCP

../../../Cedric_Memoire_Code/coexCurveFS/coexCurveFS \
	../input_files/input_coexCurveFS.dat
mv log.dat logs/log_coexCurveFS.dat

cd ..

### get FCC coexCurves

cp -r ../WHDF_rc12_diagram/coexCurves/* ./coexCurves_FCC/

### plot diagram

cd diagram
./plot_diagram.sh
cd ..