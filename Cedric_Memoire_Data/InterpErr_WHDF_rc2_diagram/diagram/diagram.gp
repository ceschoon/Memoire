#----------------------------------------
set terminal epslatex standalone color size 9cm,7cm
set output 'diagram.tex'

set title "WHDF Phase diagram, $r_c=2$" font ",20"
set xlabel '$\rho$'
set ylabel '$kT$'

set style data lines
set grid
set key bottom left

set xrange[0:1.0]

set label 'V' at 0.02,1.04
set label 'L' at 0.62,1.04
set label 'S' at 0.92,1.04

plot '../coexCurves/coexCurveVL.dat' using 2:1 linecolor 1 title 'DFT model', \
     '../coexCurves/coexCurveVL.dat' using 3:1 linecolor 1 notitle , \
     '../coexCurves/coexCurveVS.dat' using 2:1 linecolor 1 notitle , \
     '../coexCurves/coexCurveVS.dat' using 3:1 linecolor 1 notitle , \
     '../coexCurves/coexCurveLS.dat' using 2:1 linecolor 1 notitle , \
     '../coexCurves/coexCurveLS.dat' using 3:1 linecolor 1 notitle 