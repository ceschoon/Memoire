set terminal epslatex standalone color size 18cm,5cm
set output 'multiplot.tex'

set multiplot layout 1,3


##################################

set key off
unset xtics
unset ytics
set xzeroaxis linetype 1 linecolor 8
set yzeroaxis linetype 1 linecolor 8
set style data lines


######### Plot potentials #########

set xlabel '\Large$r$'
set ylabel '\Large$v(r)$'
#set xrange [-0.2:2.2]
set xrange [0.8:2.2]
set yrange [-1.2:1.2]

plot "../data/potentials.dat" using 1:3 with lines linecolor 1 linewidth 3 notitle ,\
     "../data/potentials.dat" using 1:4 with lines linecolor 2 linewidth 3 notitle


######### Plot diagram 1 #########

set xlabel '$\rho$'
set ylabel '$kT$'
set xrange[0:1.0]
set yrange [0.4:1.4]

set label 1 'V' at 0.02,1.04
set label 2 'L' at 0.62,1.04
set label 3 'S' at 0.92,1.04

plot '../data/coexCurveVL_WHDF_rc2.dat' using 2:1 linecolor 1 linewidth 3 notitle , \
     '../data/coexCurveVL_WHDF_rc2.dat' using 3:1 linecolor 1 linewidth 3 notitle , \
     '../data/coexCurveVS_WHDF_rc2.dat' using 2:1 linecolor 1 linewidth 3 notitle , \
     '../data/coexCurveVS_WHDF_rc2.dat' using 3:1 linecolor 1 linewidth 3 notitle , \
     '../data/coexCurveLS_WHDF_rc2.dat' using 2:1 linecolor 1 linewidth 3 notitle , \
     '../data/coexCurveLS_WHDF_rc2.dat' using 3:1 linecolor 1 linewidth 3 notitle 

unset label 1
unset label 2
unset label 3

######### Plot diagram 2 #########

set xlabel '$\rho$'
set ylabel '$kT$'
set xrange[0:1.3]
set yrange [0.3:1.0]

set label 4 'F' at 0.28,0.75
set label 5 'S' at 1.22,0.75

plot '../data/coexCurveFS_WHDF_rc12.dat' using 2:1 linecolor 2 linewidth 3 notitle , \
     '../data/coexCurveFS_WHDF_rc12.dat' using 3:1 linecolor 2 linewidth 3 notitle , \
#     '../data/coexCurveVL_WHDF_rc12.dat' using 2:1 linecolor 2 linewidth 3 notitle , \
#     '../data/coexCurveVL_WHDF_rc12.dat' using 3:1 linecolor 2 linewidth 3 notitle 

unset label 4
unset label 5

##################################

unset multiplot

