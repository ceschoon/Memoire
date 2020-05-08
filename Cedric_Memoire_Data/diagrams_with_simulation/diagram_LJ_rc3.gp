#----------------------------------------
set terminal epslatex standalone color size 9cm,7cm
set output 'diagram_LJ_rc3.tex'

set title "LJ Phase diagram, $r_c=3$" font ",20"
set xlabel '$\rho$'
set ylabel '$kT$'

set style data lines
set grid
set key bottom left

set xrange[0:1.1]
set yrange[0:1.6]

set label 'V' at 0.00,1.1
set label 'L' at 0.72,1.1
set label 'S' at 1.02,0.7

plot 'LJ_rc3_data_dft/coexCurveVL.dat' using 2:1 linecolor 1 title 'DFT model', \
     'LJ_rc3_data_dft/coexCurveVL.dat' using 3:1 linecolor 1 notitle , \
     'LJ_rc3_data_dft/coexCurveVS.dat' using 2:1 linecolor 1 notitle , \
     'LJ_rc3_data_dft/coexCurveVS.dat' using 3:1 linecolor 1 notitle , \
     'LJ_rc3_data_dft/coexCurveLS.dat' using 2:1 linecolor 1 notitle , \
     'LJ_rc3_data_dft/coexCurveLS.dat' using 3:1 linecolor 1 notitle , \
     'LJ_rc3_data_sim/Schultz_Coex_VL_V.dat' linecolor 2 title 'simulation', \
     'LJ_rc3_data_sim/Schultz_Coex_VL_L.dat' linecolor 2 notitle , \
     'LJ_rc3_data_sim/Schultz_Coex_VS_V.dat' linecolor 2 notitle , \
     'LJ_rc3_data_sim/Schultz_Coex_VS_S.dat' linecolor 2 notitle , \
     'LJ_rc3_data_sim/Schultz_Coex_LS_L.dat' linecolor 2 notitle , \
     'LJ_rc3_data_sim/Schultz_Coex_LS_S.dat' linecolor 2 notitle 

