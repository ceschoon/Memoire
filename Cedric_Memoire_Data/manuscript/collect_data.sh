#! /bin/bash

mkdir -p data

# collect coexistence curves data (includes vacancy data)

cp ../LJ_rc3_diagram/coexCurves/coexCurveLS.dat     data/coexCurveLS_LJ_rc3.dat
cp ../LJ_rc3_diagram/coexCurves/coexCurveVS.dat     data/coexCurveVS_LJ_rc3.dat
cp ../LJ_rc4_diagram/coexCurves/coexCurveLS.dat     data/coexCurveLS_LJ_rc4.dat
cp ../LJ_rc4_diagram/coexCurves/coexCurveVS.dat     data/coexCurveVS_LJ_rc4.dat
cp ../WHDF_rc2_diagram/coexCurves/coexCurveLS.dat   data/coexCurveLS_WHDF_rc2.dat
cp ../WHDF_rc2_diagram/coexCurves/coexCurveVS.dat   data/coexCurveVS_WHDF_rc2.dat
cp ../WHDF_rc12_diagram_Analytic/coexCurves/coexCurveFS.dat  data/coexCurveFS_WHDF_rc12.dat

# combine LJ rc=4 results from runs at higher temperature 

cd data
cat coexCurveLS_LJ_rc4.dat ../../vacancies/data_dft/coexCurveLS_LJ_rc4_highT.dat > coexCurveLS_LJ_rc4_long.dat
cd ..

# collect LJ simulation data for vacancy concentration

touch data/simulation_vacancies_LJ.dat
echo "# Note: The vacancy concentration is 'phi'" > data/simulation_vacancies_LJ.dat
echo "# " >> data/simulation_vacancies_LJ.dat
cat ../vacancies/data_sim/out_purohit2018_coex.dat  >> data/simulation_vacancies_LJ.dat

# collect older DFT results for vacancy concentration

cp ../vacancies/data_other_dft/mcrae1990.dat   data/mcrae1990.dat
cp ../vacancies/data_other_dft/singh2007.dat   data/singh2007.dat

# collect FCC-HCP free energy comparison data

cp ../WHDF_rc2_comparison_FCC_HCP/fixedkTMu/plot/data.dat    data/FCCvsHCP_WHDF_rc2.dat
cp ../WHDF_rc12_comparison_FCC_HCP_Analytic/fixedkTMu/plot/data.dat   data/FCCvsHCP_WHDF_rc12.dat

