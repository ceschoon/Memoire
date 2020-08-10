#----------------------------------------
set terminal epslatex standalone color size 9cm,7cm
set output 'diagram.tex'

set title "FCC vs HCP (WHDF $r_c=1.2$)" font ",20"
set xlabel '$\rho$'
set ylabel '$kT$'

set style data lines
set grid
set key top left

set xrange[0:1.3]

set label 'F' at 0.28,0.75
set label 'S' at 1.22,0.75

plot '../coexCurves_FCC/coexCurveVL.dat' using 2:1 linecolor 1 notitle, \
     '../coexCurves_FCC/coexCurveVL.dat' using 3:1 linecolor 1 notitle , \
     '../coexCurves_FCC/coexCurveFS.dat' using 2:1 linecolor 1 title 'FCC' , \
     '../coexCurves_FCC/coexCurveFS.dat' using 3:1 linecolor 1 notitle , \
     '../coexCurves_HCP/coexCurveFS.dat' using 2:1 with points pointtype 7 linecolor 2 title 'HCP' , \
     '../coexCurves_HCP/coexCurveFS.dat' using 3:1 with points pointtype 7 linecolor 2 notitle , \