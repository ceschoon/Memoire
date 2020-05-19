#----------------------------------------
set terminal epslatex standalone color size 6cm,5cm
set output 'diagram_WHDF_rc12.tex'

#set title "WHDF Phase diagram, $r_c=1.2$" font ",20"
set xlabel '$\rho$'
set ylabel '$kT$'

set style data lines
set key off

set xrange[0:1.3]
unset xtics
unset ytics

set label 'F' at 0.28,0.75
set label 'S' at 1.22,0.75

plot '../data/coexCurveVL_WHDF_rc12.dat' using 2:1 linecolor 2 notitle , \
     '../data/coexCurveVL_WHDF_rc12.dat' using 3:1 linecolor 2 notitle , \
     '../data/coexCurveFS_WHDF_rc12.dat' using 2:1 linecolor 2 notitle , \
     '../data/coexCurveFS_WHDF_rc12.dat' using 3:1 linecolor 2 notitle 