#----------------------------------------
set terminal epslatex standalone color size 9cm,7cm
set output 'diagram.tex'

set title "FCC vs HCP (WHDF $r_c=2$)" font ",20"
set xlabel '$\rho$'
set ylabel '$kT$'

set style data lines
set grid
set key at 0.5,0.68 #bottom left

set xrange[0:1.0]

set label 'V' at 0.02,1.04
set label 'L' at 0.62,1.04
set label 'S' at 0.92,1.04

plot '../coexCurves_FCC/coexCurveVL.dat' using 2:1 linecolor 1 notitle, \
     '../coexCurves_FCC/coexCurveVL.dat' using 3:1 linecolor 1 notitle , \
     '../coexCurves_FCC/coexCurveVS.dat' using 2:1 linecolor 1 title 'FCC' , \
     '../coexCurves_FCC/coexCurveVS.dat' using 3:1 linecolor 1 notitle , \
     '../coexCurves_FCC/coexCurveLS.dat' using 2:1 linecolor 1 notitle , \
     '../coexCurves_FCC/coexCurveLS.dat' using 3:1 linecolor 1 notitle , \
     '../coexCurves_HCP/coexCurveVS.dat' using 2:1 with points pointtype 7 linecolor 2 title 'HCP' , \
     '../coexCurves_HCP/coexCurveVS.dat' using 3:1 with points pointtype 7 linecolor 2 notitle , \
     '../coexCurves_HCP/coexCurveLS.dat' using 2:1 with points pointtype 7 linecolor 2 notitle , \
     '../coexCurves_HCP/coexCurveLS.dat' using 3:1 with points pointtype 7 linecolor 2 notitle 