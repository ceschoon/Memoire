#! /bin/bash

mkdir -p data

cp ../TEST/potentials.dat data/

cp ../../../Cedric_Memoire_Data/WHDF_rc2_diagram/coexCurves/coexCurveVL.dat \
	data/coexCurveVL_WHDF_rc2.dat
cp ../../../Cedric_Memoire_Data/WHDF_rc2_diagram/coexCurves/coexCurveVS.dat \
	data/coexCurveVS_WHDF_rc2.dat
cp ../../../Cedric_Memoire_Data/WHDF_rc2_diagram/coexCurves/coexCurveLS.dat \
	data/coexCurveLS_WHDF_rc2.dat

cp ../../../Cedric_Memoire_Data/WHDF_rc12_diagram/coexCurves/coexCurveVL.dat \
	data/coexCurveVL_WHDF_rc12.dat
cp ../../../Cedric_Memoire_Data/WHDF_rc12_diagram/coexCurves/coexCurveFS.dat \
	data/coexCurveFS_WHDF_rc12.dat

