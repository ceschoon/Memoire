#----------------------------------------
set terminal epslatex standalone color size 9cm,7cm
set output 'diagram.tex'

set title "WHDF Phase diagram, $r_c=1.2$" font ",20"
set xlabel '$\rho$'
set ylabel '$kT$'

set style data lines
set grid
set key top left

set label 'F' at 0.3,0.8
set label 'S' at 1.1,0.8

plot '../coexCurves/coexCurveVL.dat' using 2:1 linecolor 1 title 'DFT model', \
     '../coexCurves/coexCurveVL.dat' using 3:1 linecolor 1 notitle , \
     '../coexCurves/coexCurveFS.dat' using 2:1 linecolor 1 notitle , \
     '../coexCurves/coexCurveFS.dat' using 3:1 linecolor 1 notitle 