#! /bin/bash

outFile=summary.txt

rm $outFile
touch $outFile

echo "" >> $outFile
echo "# Errors in HCP interpolation (WHDF rc=2.0)" >> $outFile
cat ../InterpErr_WHDF_rc2_comparison_FCC_HCP_Analytic/fixedkTMu/dataHCP/kTMuSolid* | grep err \
	> HCP_rc2.temp
cat HCP_rc2.temp >> $outFile
echo $(./compute_mean HCP_rc2.temp) >> $outFile
echo "" >> $outFile

echo "" >> $outFile
echo "# Errors in HCP interpolation (WHDF rc=1.2)" >> $outFile
cat ../InterpErr_WHDF_rc12_comparison_FCC_HCP_Analytic/fixedkTMu/dataHCP/kTMuSolid* | grep err \
	> HCP_rc12.temp
cat HCP_rc12.temp >> $outFile
echo $(./compute_mean HCP_rc12.temp) >> $outFile
echo "" >> $outFile

echo "" >> $outFile
echo "# Errors in min wrt Ngrid (WHDF rc=2.0 coex VS) (print a few)" >> $outFile
cat ../InterpErr_WHDF_rc2_diagram/dataVS/kTMuSolid* | grep err > coexVS_rc2.temp
head coexVS_rc2.temp >> $outFile
echo $(./compute_mean coexVS_rc2.temp) >> $outFile
echo "" >> $outFile

echo "" >> $outFile
echo "# Errors in min wrt Ngrid (WHDF rc=2.0 coex LS) (print a few)" >> $outFile
cat ../InterpErr_WHDF_rc2_diagram/dataLS/kTMuSolid* | grep err > coexLS_rc2.temp
head coexLS_rc2.temp >> $outFile
echo $(./compute_mean coexLS_rc2.temp) >> $outFile
echo "" >> $outFile

rm *.temp
cat $outFile