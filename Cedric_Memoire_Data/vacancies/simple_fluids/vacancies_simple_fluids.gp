
set term epslatex standalone color size 9cm,7cm
set output 'vacancies_simple_fluids.tex'

set title "Vacancies along coexistence curves"
set xlabel "$kT$"
set ylabel "$c_{vac} \\times 10^{-4}$"
set style data lines
set grid linecolor rgb '#B0B0B0'
#set key at 3.7,1.6 
set key top right

set xrange[0.4:1.8]
set yrange[0:10]

plot '../data_dft/coexCurveVS_LJ_rc3.dat' using 1:($7*1e4) with lines linecolor 1 notitle, \
     '../data_dft/coexCurveLS_LJ_rc3.dat' using 1:($7*1e4) with lines linecolor 1 title 'LJ $r_c=3$', \
     '../data_dft/coexCurveVS_LJ_rc4.dat' using 1:(log($7)/log(10)) with lines linecolor 8 notitle, \
     '../data_dft/coexCurveLS_LJ_rc4.dat' using 1:(log($7)/log(10)) with lines linecolor 8 title 'LJ $r_c=4$', \
     '../data_dft/coexCurveVS_WHDF_rc2.dat' using 1:($7*1e4) with lines linecolor 2 notitle, \
     '../data_dft/coexCurveLS_WHDF_rc2.dat' using 1:($7*1e4) with lines linecolor 2 title 'WHDF $r_c=2$', \
