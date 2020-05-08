#----------------------------------------
set terminal epslatex standalone color size 9cm,7cm
set output 'diagram_WHDF_rc2.tex'

set title "WHDF Phase diagram, $r_c=2$" font ",20"
set xlabel '$\rho$'
set ylabel '$kT$'

set style data lines
set grid
set key bottom left

set xrange[0:1.0]
set yrange[0:1.4]

set label 'V' at 0.02,1.1
set label 'L' at 0.62,1.1
set label 'S' at 0.95,0.5

plot 'WHDF_rc2_data_dft/coexCurveVL.dat' using 2:1 linecolor 1 title 'DFT model', \
     'WHDF_rc2_data_dft/coexCurveVL.dat' using 3:1 linecolor 1 notitle , \
     'WHDF_rc2_data_dft/coexCurveVS.dat' using 2:1 linecolor 1 notitle , \
     'WHDF_rc2_data_dft/coexCurveVS.dat' using 3:1 linecolor 1 notitle , \
     'WHDF_rc2_data_dft/coexCurveLS.dat' using 2:1 linecolor 1 notitle , \
     'WHDF_rc2_data_dft/coexCurveLS.dat' using 3:1 linecolor 1 notitle , \
     'WHDF_rc2_data_sim/coexCurveVL_Frenkel.dat' using 2:1 linecolor 2 title 'simulation', \
     'WHDF_rc2_data_sim/coexCurveVL_Frenkel.dat' using 3:1 linecolor 2 notitle , \
     'WHDF_rc2_data_sim/coexCurveVS_Frenkel.dat' using 2:1 linecolor 2 notitle , \
     'WHDF_rc2_data_sim/coexCurveVS_Frenkel.dat' using 3:1 linecolor 2 notitle , \
     'WHDF_rc2_data_sim/coexCurveLS_Frenkel.dat' using 2:1 linecolor 2 notitle , \
     'WHDF_rc2_data_sim/coexCurveLS_Frenkel.dat' using 3:1 linecolor 2 notitle 

