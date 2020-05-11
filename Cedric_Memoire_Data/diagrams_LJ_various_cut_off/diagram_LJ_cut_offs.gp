#----------------------------------------
set terminal epslatex standalone color size 9cm,7cm
set output 'diagram_LJ_cut_offs.tex'

set title "Effect of the cut-off (LJ)" font ",20"
set xlabel '$\rho$'
set ylabel '$kT$'

set style data lines
set grid
set key bottom left

set xrange[0:1.1]
set yrange[0:1.6]

set label 'V' at 0.02,1.3
set label 'L' at 0.68,1.3
set label 'S' at 1.03,0.9

plot 'LJ_rc3_data_dft/coexCurveVL.dat' using 2:1 linecolor 1 title '$r_c=3$', \
     'LJ_rc3_data_dft/coexCurveVL.dat' using 3:1 linecolor 1 notitle , \
     'LJ_rc3_data_dft/coexCurveVS.dat' using 2:1 linecolor 1 notitle , \
     'LJ_rc3_data_dft/coexCurveVS.dat' using 3:1 linecolor 1 notitle , \
     'LJ_rc3_data_dft/coexCurveLS.dat' using 2:1 linecolor 1 notitle , \
     'LJ_rc3_data_dft/coexCurveLS.dat' using 3:1 linecolor 1 notitle , \
     'LJ_rc4_data_dft/coexCurveVL.dat' using 2:1 linecolor 2 title '$r_c=4$', \
     'LJ_rc4_data_dft/coexCurveVL.dat' using 3:1 linecolor 2 notitle , \
     'LJ_rc4_data_dft/coexCurveVS.dat' using 2:1 linecolor 2 notitle , \
     'LJ_rc4_data_dft/coexCurveVS.dat' using 3:1 linecolor 2 notitle , \
     'LJ_rc4_data_dft/coexCurveLS.dat' using 2:1 linecolor 2 notitle , \
     'LJ_rc4_data_dft/coexCurveLS.dat' using 3:1 linecolor 2 notitle , \
     'LJ_rc10_data_dft/coexCurveVL.dat' using 2:1 linecolor 8 title '$r_c=10$', \
     'LJ_rc10_data_dft/coexCurveVL.dat' using 3:1 linecolor 8 notitle , \
     'LJ_rc10_data_dft/coexCurveVS.dat' using 2:1 linecolor 8 notitle , \
     'LJ_rc10_data_dft/coexCurveVS.dat' using 3:1 linecolor 8 notitle , \
     'LJ_rc10_data_dft/coexCurveLS.dat' using 2:1 linecolor 8 notitle , \
     'LJ_rc10_data_dft/coexCurveLS.dat' using 3:1 linecolor 8 notitle 
