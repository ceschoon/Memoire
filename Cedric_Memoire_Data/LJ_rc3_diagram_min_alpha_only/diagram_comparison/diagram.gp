#----------------------------------------
set terminal epslatex standalone color size 9cm,7cm
set output 'diagram.tex'

set title "No minimisation over vacancies (LJ $r_c=3$)" font ",20"
set xlabel '$\rho$'
set ylabel '$kT$'

set style data lines
set grid
set key bottom left

set xrange[0:1.1]

set label 'V' at 0.01,1.1
set label 'L' at 0.62,1.1
set label 'S' at 1.02,1.1

plot '../coexCurves/coexCurveVL.dat' using 2:1 linecolor 2 title 'without $c_{vac}$', \
     '../coexCurves/coexCurveVL.dat' using 3:1 linecolor 2 notitle , \
     '../coexCurves/coexCurveVS.dat' using 2:1 linecolor 2 notitle , \
     '../coexCurves/coexCurveVS.dat' using 3:1 linecolor 2 notitle , \
     '../coexCurves/coexCurveLS.dat' using 2:1 linecolor 2 notitle , \
     '../coexCurves/coexCurveLS.dat' using 3:1 linecolor 2 notitle , \
     '../coexCurves_with_Cvac/coexCurveVL.dat' using 2:1 linecolor 1 title 'with $c_{vac}$', \
     '../coexCurves_with_Cvac/coexCurveVL.dat' using 3:1 linecolor 1 notitle , \
     '../coexCurves_with_Cvac/coexCurveVS.dat' using 2:1 linecolor 1 notitle , \
     '../coexCurves_with_Cvac/coexCurveVS.dat' using 3:1 linecolor 1 notitle , \
     '../coexCurves_with_Cvac/coexCurveLS.dat' using 2:1 linecolor 1 notitle , \
     '../coexCurves_with_Cvac/coexCurveLS.dat' using 3:1 linecolor 1 notitle , \
