#! /bin/bash

cp ../LJ_rc3_diagram/coexCurves/coexCurveLS.dat     data_dft/coexCurveLS_LJ_rc3.dat
cp ../LJ_rc3_diagram/coexCurves/coexCurveVS.dat     data_dft/coexCurveVS_LJ_rc3.dat
cp ../LJ_rc4_diagram/coexCurves/coexCurveLS.dat     data_dft/coexCurveLS_LJ_rc4.dat
cp ../LJ_rc4_diagram/coexCurves/coexCurveVS.dat     data_dft/coexCurveVS_LJ_rc4.dat
cp ../WHDF_rc2_diagram/coexCurves/coexCurveLS.dat   data_dft/coexCurveLS_WHDF_rc2.dat
cp ../WHDF_rc2_diagram/coexCurves/coexCurveVS.dat   data_dft/coexCurveVS_WHDF_rc2.dat
cp ../WHDF_rc12_diagram/coexCurves/coexCurveFS.dat  data_dft/coexCurveFS_WHDF_rc12.dat

# combine LJ rc=4 results from runs at higher temperature 
cd data_dft
cat coexCurveLS_LJ_rc4.dat coexCurveLS_LJ_rc4_highT.dat > coexCurveLS_LJ_rc4_long.dat
cd ..

