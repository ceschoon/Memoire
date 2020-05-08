#----------------------------------------
set terminal epslatex standalone color size 9cm,7cm
set output 'diagram_WHDF_rc12.tex'

set title "WHDF Phase diagram, $r_c=1.2$" font ",20"
set xlabel '$\rho$'
set ylabel '$kT$'

set style data lines
set grid
set key bottom left

set xrange[0:1.3]
set yrange[0:1.0]

set label 'F' at 0.28,0.7
set label 'S' at 1.22,0.7

plot 'WHDF_rc12_data_dft/coexCurveVL.dat' using 2:1 linecolor 1 title 'DFT model', \
     'WHDF_rc12_data_dft/coexCurveVL.dat' using 3:1 linecolor 1 notitle , \
     'WHDF_rc12_data_dft/coexCurveFS.dat' using 2:1 linecolor 1 notitle , \
     'WHDF_rc12_data_dft/coexCurveFS.dat' using 3:1 linecolor 1 notitle , \
     'WHDF_rc12_data_sim/coexCurveFS_Frenkel.dat' using 2:1 linecolor 2 title 'simulation', \
     'WHDF_rc12_data_sim/coexCurveFS_Frenkel.dat' using 3:1 linecolor 2 notitle 

