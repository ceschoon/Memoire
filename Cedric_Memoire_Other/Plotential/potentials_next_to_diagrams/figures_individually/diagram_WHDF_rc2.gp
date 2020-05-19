#----------------------------------------
set terminal epslatex standalone color size 6cm,5cm
set output 'diagram_WHDF_rc2.tex'

#set title "WHDF Phase diagram, $r_c=2$" font ",20"
set xlabel '$\rho$'
set ylabel '$kT$'

set style data lines
set key off

set xrange[0:1.0]
unset xtics
unset ytics

set label 'V' at 0.02,1.04
set label 'L' at 0.62,1.04
set label 'S' at 0.92,1.04

plot '../data/coexCurveVL_WHDF_rc2.dat' using 2:1 linecolor 1 notitle , \
     '../data/coexCurveVL_WHDF_rc2.dat' using 3:1 linecolor 1 notitle , \
     '../data/coexCurveVS_WHDF_rc2.dat' using 2:1 linecolor 1 notitle , \
     '../data/coexCurveVS_WHDF_rc2.dat' using 3:1 linecolor 1 notitle , \
     '../data/coexCurveLS_WHDF_rc2.dat' using 2:1 linecolor 1 notitle , \
     '../data/coexCurveLS_WHDF_rc2.dat' using 3:1 linecolor 1 notitle 