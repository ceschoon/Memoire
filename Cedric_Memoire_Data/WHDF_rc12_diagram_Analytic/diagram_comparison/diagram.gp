#----------------------------------------
set terminal epslatex standalone color size 9cm,7cm
set output 'diagram.tex'

set title "Analytic vs Numeric Weights (WHDF $r_c=1.2$) $\\quad$"
set xlabel '$\rho$'
set ylabel '$kT$'

set style data lines
set grid
set key top left

set xrange[0:1.3]

set label 'F' at 0.28,0.75
set label 'S' at 1.22,0.75

plot '../coexCurves/coexCurveVL.dat' using 2:1 linecolor 1 title 'Analytic', \
     '../coexCurves/coexCurveVL.dat' using 3:1 linecolor 1 notitle , \
     '../coexCurves/coexCurveFS.dat' using 2:1 linecolor 1 notitle , \
     '../coexCurves/coexCurveFS.dat' using 3:1 linecolor 1 notitle , \
     '../coexCurves_Num/coexCurveVL.dat' using 2:1 linecolor 2 title 'Numeric', \
     '../coexCurves_Num/coexCurveVL.dat' using 3:1 linecolor 2 notitle , \
     '../coexCurves_Num/coexCurveFS.dat' using 2:1 linecolor 2 notitle , \
     '../coexCurves_Num/coexCurveFS.dat' using 3:1 linecolor 2 notitle , \

