#! /bin/bash

mkdir -p LJ_rc3_data_dft
mkdir -p WHDF_rc2_data_dft
mkdir -p WHDF_rc12_data_dft

mkdir -p LJ_rc3_data_sim
mkdir -p WHDF_rc2_data_sim
mkdir -p WHDF_rc12_data_sim

cp ../LJ_rc3_diagram/coexCurves/coexCurve*.dat LJ_rc3_data_dft/
cp ../WHDF_rc2_diagram/coexCurves/coexCurve*.dat WHDF_rc2_data_dft/
cp ../WHDF_rc12_diagram/coexCurves/coexCurve*.dat WHDF_rc12_data_dft/

cp ../../Cedric_Memoire_Other/Coex_Frenkel/coexCurveVL_Frenkel.dat WHDF_rc2_data_sim/
cp ../../Cedric_Memoire_Other/Coex_Frenkel/coexCurveVS_Frenkel.dat WHDF_rc2_data_sim/
cp ../../Cedric_Memoire_Other/Coex_Frenkel/coexCurveLS_Frenkel.dat WHDF_rc2_data_sim/
cp ../../Cedric_Memoire_Other/Coex_Frenkel/coexCurveFS_Frenkel.dat WHDF_rc12_data_sim/

cp ../../Cedric_Memoire_Other/Coex_Schultz_wpd/Schultz_Coex_*.dat LJ_rc3_data_sim/
