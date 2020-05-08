#----------------------------------------
set terminal epslatex standalone color size 9cm,7cm
set output 'diagram_example.tex'

set title "Phase diagram: example" font ",20"
set xlabel '$\rho$'
set ylabel '$kT$'

set style data lines
set grid
set key bottom left

set xrange[0:1.1]
set yrange[0.5:1.5]

set label 'V' at 0.01,1.1
set label 'L' at 0.62,1.1
set label 'S' at 1.02,1.1
set label 'F' at 0.28,1.44
set label 'U' at 0.25,0.9
set label 'U' at 0.80,0.7

set label at 0.244, 1.287 "" point pointtype 7 linecolor 4 front
set label at 0.039, 1.000 "" point pointtype 7 linecolor 7
set label at 0.561, 1.000 "" point pointtype 7 linecolor 7
set label at 0.001, 0.600 "" point pointtype 7 linecolor 1
set label at 0.962, 0.600 "" point pointtype 7 linecolor 1
set label at 0.774, 1.000 "" point pointtype 7 linecolor 2
set label at 0.922, 1.000 "" point pointtype 7 linecolor 2

set arrow from 0.0,0.792 to 1.1,0.792 nohead linewidth 2 dashtype 2 linecolor 4 

plot 'LJ_rc3_data_dft/coexCurveVL.dat' using 2:1 linecolor 7 notitle , \
     'LJ_rc3_data_dft/coexCurveVL.dat' using 3:1 linecolor 7 notitle , \
     'LJ_rc3_data_dft/coexCurveVS.dat' using 2:1 linecolor 1 notitle , \
     'LJ_rc3_data_dft/coexCurveVS.dat' using 3:1 linecolor 1 notitle , \
     'LJ_rc3_data_dft/coexCurveLS.dat' using 2:1 linecolor 2 notitle , \
     'LJ_rc3_data_dft/coexCurveLS.dat' using 3:1 linecolor 2 notitle 



