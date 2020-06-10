#----------------------------------------
set terminal svg size 600,500
set output 'diagram_WHDF_rc12.svg'

#set title "WHDF Phase diagram, $r_c=1.2$" font ",30"
set xlabel 'density' font ",30" #offset 0,-1
set ylabel 'temperature' font ",30" #offset -1,0

set style data lines
#set grid
#set key bottom left

set xrange[0.0:1.5]
set yrange[0.3:1.0]

#set xtics font ",30"
#set ytics font ",30"
unset xtics
unset ytics

set label 'Fluid' at 0.28,0.8 font ",30" front
set label 'Solid' at 1.22,0.7 font ",30" front
set label 'Unstable' at 0.6,0.5 font ",30" front

plot 'WHDF_rc12_data_dft/coexCurveVL.dat' using 2:1 lc 2 lw 2 dt '-_' notitle, \
     'WHDF_rc12_data_dft/coexCurveVL.dat' using 3:1 lc 2 lw 2 dt '-_' notitle , \
     'WHDF_rc12_data_dft/coexCurveFS.dat' using 2:1 lc 2 lw 2 notitle , \
     'WHDF_rc12_data_dft/coexCurveFS.dat' using 3:1 lc 2 lw 2 notitle , \

