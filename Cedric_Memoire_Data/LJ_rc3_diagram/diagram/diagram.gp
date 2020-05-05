#----------------------------------------
set terminal epslatex standalone color size 9cm,7cm
set output 'diagram.tex'

set title "LJ Phase diagram, $r_c=3$" font ",20"
set xlabel '$\rho$'
set ylabel '$kT$'

set style data lines
set grid
set key bottom left

set label 'V' at 0.03,1.2
set label 'L' at 0.60,1.2
set label 'S' at 1.00,1.2

plot '../coexCurves/coexCurveVL.dat' using 2:1 linecolor 1 title 'DFT model', \
     '../coexCurves/coexCurveVL.dat' using 3:1 linecolor 1 notitle , \
     '../coexCurves/coexCurveVS.dat' using 2:1 linecolor 1 notitle , \
     '../coexCurves/coexCurveVS.dat' using 3:1 linecolor 1 notitle , \
     '../coexCurves/coexCurveLS.dat' using 2:1 linecolor 1 notitle , \
     '../coexCurves/coexCurveLS.dat' using 3:1 linecolor 1 notitle 