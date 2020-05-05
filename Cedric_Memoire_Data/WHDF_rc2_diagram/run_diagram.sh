#! /bin/bash

mkdir -p coexCurves
mkdir -p coexCurves/logs
mkdir -p critical_triple_points

### generate coexistence curves

cd coexCurves

../../../Cedric_Memoire_Code/coexCurveVL/coexCurveVL \
	../input_files/input_coexCurveVL.dat
mv log.dat logs/log_coexCurveVL.dat

../../../Cedric_Memoire_Code/coexCurveFS/coexCurveFS \
	../input_files/input_coexCurveVS.dat
mv log.dat logs/log_coexCurveVS.dat

../../../Cedric_Memoire_Code/coexCurveFS/coexCurveFS \
	../input_files/input_coexCurveLS.dat
mv log.dat logs/log_coexCurveLS.dat

cd ..

### draw phase diagram

cd diagram
./plot_diagram.sh
cd ..

### compute critical point

cd critical_triple_points
../../../Cedric_Memoire_Code/criticalPoint/criticalPoint \
	../input_files/input_criticalPoint.dat
mv log.dat log_criticalPoint.dat
cd ..

### fits of fluid-solid coexistence curves

cp coexCurves/coexCurve*.dat fits/data/

cd fits
../../../Cedric_Memoire_Code/polyFit/polyFit input_VS_mu.dat
../../../Cedric_Memoire_Code/polyFit/polyFit input_VS_omega.dat
../../../Cedric_Memoire_Code/polyFit/polyFit input_VS_rhoV.dat
../../../Cedric_Memoire_Code/polyFit/polyFit input_VS_rhoS.dat

../../../Cedric_Memoire_Code/polyFit/polyFit input_LS_mu.dat
../../../Cedric_Memoire_Code/polyFit/polyFit input_LS_omega.dat
../../../Cedric_Memoire_Code/polyFit/polyFit input_LS_rhoL.dat
../../../Cedric_Memoire_Code/polyFit/polyFit input_LS_rhoS.dat
cd ..

### compute triple point

cd critical_triple_points
../../../Cedric_Memoire_Code/triplePoint/triplePoint \
	../input_files/input_triplePoint.dat
mv log.dat log_triplePoint.dat
cd ..
