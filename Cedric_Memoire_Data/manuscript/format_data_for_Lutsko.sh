#! /bin/bash

mkdir -p data_formatted





###########################
# vacancies in gaussian DFT
###########################

# vacancies for LJ rc=3

./code_formatting/extract_one_column data/coexCurveVS_LJ_rc3.dat 0 \
	> data_formatted/vacancies_LJ_rc3_VS_x.temp

./code_formatting/extract_one_column data/coexCurveVS_LJ_rc3.dat 6 \
	> data_formatted/vacancies_LJ_rc3_VS_y.temp

./code_formatting/extract_one_column data/coexCurveLS_LJ_rc3.dat 0 \
	> data_formatted/vacancies_LJ_rc3_LS_x.temp

./code_formatting/extract_one_column data/coexCurveLS_LJ_rc3.dat 6 \
	> data_formatted/vacancies_LJ_rc3_LS_y.temp

# vacancies for LJ rc=4

./code_formatting/extract_one_column data/coexCurveVS_LJ_rc4.dat 0 \
	> data_formatted/vacancies_LJ_rc4_VS_x.temp

./code_formatting/extract_one_column data/coexCurveVS_LJ_rc4.dat 6 \
	> data_formatted/vacancies_LJ_rc4_VS_y.temp

./code_formatting/extract_one_column data/coexCurveLS_LJ_rc4.dat 0 \
	> data_formatted/vacancies_LJ_rc4_LS_x.temp

./code_formatting/extract_one_column data/coexCurveLS_LJ_rc4.dat 6 \
	> data_formatted/vacancies_LJ_rc4_LS_y.temp

# vacancies for WHDF rc=2

./code_formatting/extract_one_column data/coexCurveVS_WHDF_rc2.dat 0 \
	> data_formatted/vacancies_WHDF_rc2_VS_x.temp

./code_formatting/extract_one_column data/coexCurveVS_WHDF_rc2.dat 6 \
	> data_formatted/vacancies_WHDF_rc2_VS_y.temp

./code_formatting/extract_one_column data/coexCurveLS_WHDF_rc2.dat 0 \
	> data_formatted/vacancies_WHDF_rc2_LS_x.temp

./code_formatting/extract_one_column data/coexCurveLS_WHDF_rc2.dat 6 \
	> data_formatted/vacancies_WHDF_rc2_LS_y.temp

# vacancies for WHDF rc=1.2

./code_formatting/extract_one_column data/coexCurveFS_WHDF_rc12.dat 0 \
	> data_formatted/vacancies_WHDF_rc12_FS_x.temp

./code_formatting/extract_one_column data/coexCurveFS_WHDF_rc12.dat 6 \
	> data_formatted/vacancies_WHDF_rc12_FS_y.temp

# comment line

echo "# LJ_rc3_VS_kT	LJ_rc3_VS_c   	LJ_rc3_LS_kT  	LJ_rc3_LS_c   	"\
"LJ_rc4_VS_kT  	LJ_rc4_VS_c   	LJ_rc4_LS_kT  	LJ_rc4_LS_c   	"\
"WHDF_rc2_VS_kT	WHDF_rc2_VS_c 	WHDF_rc2_LS_kT	WHDF_rc2_LS_c 	"\
"WHDF_rc12_kT  	WHDF_rc12_c   	"\
	>  data_formatted/vacancies_gaussian.dat

# merge columns

paste data_formatted/vacancies_LJ_rc3_VS_x.temp \
      data_formatted/vacancies_LJ_rc3_VS_y.temp \
      data_formatted/vacancies_LJ_rc3_LS_x.temp \
      data_formatted/vacancies_LJ_rc3_LS_y.temp \
      data_formatted/vacancies_LJ_rc4_VS_x.temp \
      data_formatted/vacancies_LJ_rc4_VS_y.temp \
      data_formatted/vacancies_LJ_rc4_LS_x.temp \
      data_formatted/vacancies_LJ_rc4_LS_y.temp \
      data_formatted/vacancies_WHDF_rc2_VS_x.temp \
      data_formatted/vacancies_WHDF_rc2_VS_y.temp \
      data_formatted/vacancies_WHDF_rc2_LS_x.temp \
      data_formatted/vacancies_WHDF_rc2_LS_y.temp \
      data_formatted/vacancies_WHDF_rc12_FS_x.temp \
      data_formatted/vacancies_WHDF_rc12_FS_y.temp \
      >> data_formatted/vacancies_gaussian.dat

#./code_formatting/align_correctly data_formatted/vacancies_gaussian.dat 14

#mv data_formatted/vacancies_gaussian.dat.out \
#   data_formatted/vacancies_gaussian.dat

mv data_formatted/vacancies_gaussian.dat \
   data_formatted/Fig.vacancies_gaussian.dat






######################
# vacancies literature
######################

# vacancies for LJ rc=4 (already done above)

# vacancies from simulations

./code_formatting/extract_one_column data/simulation_vacancies_LJ.dat 3 \
	> data_formatted/vacancies_purohit2018_x.temp

./code_formatting/extract_one_column data/simulation_vacancies_LJ.dat 4 \
	> data_formatted/vacancies_purohit2018_xerr.temp

./code_formatting/extract_one_column data/simulation_vacancies_LJ.dat 5 \
	> data_formatted/vacancies_purohit2018_y.temp

./code_formatting/extract_one_column data/simulation_vacancies_LJ.dat 6 \
	> data_formatted/vacancies_purohit2018_yerr.temp

# vacancies from mcrae1990

./code_formatting/extract_one_column data/mcrae1990.dat 0 \
	> data_formatted/vacancies_mcrae1990_x.temp

./code_formatting/extract_one_column data/mcrae1990.dat 1 \
	> data_formatted/vacancies_mcrae1990_y.temp

# vacancies from singh2007

./code_formatting/extract_one_column data/singh2007.dat 0 \
	> data_formatted/vacancies_singh2007_x.temp

./code_formatting/extract_one_column data/singh2007.dat 1 \
	> data_formatted/vacancies_singh2007_y.temp

# comment line

echo "# me_LJ_rc4_kT	me_LJ_rc4_c   	sim_kT  	sim_kTerr	"\
"sim_c   	sim_cerr	mcrae_kT	mcrae_c	"\
"singh_kT	singh_c	"\
	>  data_formatted/vacancies_reflit.dat

# combine data in one file

paste data_formatted/vacancies_LJ_rc4_LS_x.temp \
      data_formatted/vacancies_LJ_rc4_LS_y.temp \
      data_formatted/vacancies_purohit2018_x.temp \
      data_formatted/vacancies_purohit2018_xerr.temp \
      data_formatted/vacancies_purohit2018_y.temp \
      data_formatted/vacancies_purohit2018_yerr.temp \
      data_formatted/vacancies_mcrae1990_x.temp \
      data_formatted/vacancies_mcrae1990_y.temp \
      data_formatted/vacancies_singh2007_x.temp \
      data_formatted/vacancies_singh2007_y.temp \
      >> data_formatted/vacancies_reflit.dat

#./code_formatting/align_correctly data_formatted/vacancies_reflit.dat 8

#mv data_formatted/vacancies_reflit.dat.out \
#   data_formatted/vacancies_reflit.dat

mv data_formatted/vacancies_reflit.dat \
   data_formatted/Fig.vacancies_reflit.dat






####################
# comparison FCC HCP
####################

# free energy HCP-FCC for WHDF rc=2

./code_formatting/extract_one_column data/FCCvsHCP_WHDF_rc2.dat 2 \
	> data_formatted/FCCvsHCP_WHDF_rc2_kT.temp

./code_formatting/extract_one_column data/FCCvsHCP_WHDF_rc2.dat 8 \
	> data_formatted/FCCvsHCP_WHDF_rc2_FCC.temp

./code_formatting/extract_one_column data/FCCvsHCP_WHDF_rc2.dat 11 \
	> data_formatted/FCCvsHCP_WHDF_rc2_HCP.temp

./code_formatting/substract_two_columns \
	data_formatted/FCCvsHCP_WHDF_rc2_HCP.temp \
	data_formatted/FCCvsHCP_WHDF_rc2_FCC.temp \
	> data_formatted/FCCvsHCP_WHDF_rc2_diff.temp

# free energy HCP-FCC for WHDF rc=1.2

./code_formatting/extract_one_column data/FCCvsHCP_WHDF_rc12.dat 2 \
	> data_formatted/FCCvsHCP_WHDF_rc12_kT.temp

./code_formatting/extract_one_column data/FCCvsHCP_WHDF_rc12.dat 8 \
	> data_formatted/FCCvsHCP_WHDF_rc12_FCC.temp

./code_formatting/extract_one_column data/FCCvsHCP_WHDF_rc12.dat 11 \
	> data_formatted/FCCvsHCP_WHDF_rc12_HCP.temp

./code_formatting/substract_two_columns \
	data_formatted/FCCvsHCP_WHDF_rc12_HCP.temp \
	data_formatted/FCCvsHCP_WHDF_rc12_FCC.temp \
	> data_formatted/FCCvsHCP_WHDF_rc12_diff.temp

# comment line

echo "# rc2_kT  	rc2_diff   	rc12_kT   	rc12_diff"\
	>  data_formatted/FCCvsHCP_WHDF.dat

# combine data in one file

paste data_formatted/FCCvsHCP_WHDF_rc2_kT.temp \
      data_formatted/FCCvsHCP_WHDF_rc2_diff.temp \
      data_formatted/FCCvsHCP_WHDF_rc12_kT.temp \
      data_formatted/FCCvsHCP_WHDF_rc12_diff.temp \
      >> data_formatted/FCCvsHCP_WHDF.dat

#./code_formatting/align_correctly data_formatted/FCCvsHCP_WHDF.dat 10

#mv data_formatted/FCCvsHCP_WHDF.dat.out \
#   data_formatted/FCCvsHCP_WHDF.dat

mv data_formatted/FCCvsHCP_WHDF.dat \
   data_formatted/Fig.FCCvsHCP_WHDF.dat







#################
# clean directory
#################

rm data_formatted/*.temp



