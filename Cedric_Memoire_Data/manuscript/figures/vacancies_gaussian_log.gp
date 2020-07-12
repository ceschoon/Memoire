
set term epslatex standalone color size 9cm,6cm
set output 'vacancies_gaussian_log.tex'

#set title "Vacancies along coexistence curves"
set xlabel "$kT$"
set ylabel "$\\log_{10}c$"

set grid linecolor rgb '#B0B0B0'
set key bottom right

set xrange[0.5:1.5]
set yrange[-5:-3]

plot '../data/coexCurveFS_WHDF_rc12.dat' using 1:(log($7)/log(10)) with lines lc 4 title 'WHDF $r_c=1.2$', \
     '../data/coexCurveVS_WHDF_rc2.dat' using 1:(log($7)/log(10)) with lines lc 2 title 'WHDF $r_c=2.0$', \
     '../data/coexCurveLS_WHDF_rc2.dat' using 1:(log($7)/log(10)) with lines lc 2 notitle, \
     '../data/coexCurveVS_LJ_rc3.dat' using 1:(log($7)/log(10)) with lines lc 7 title 'LJ $r_c=3.0$', \
     '../data/coexCurveLS_LJ_rc3.dat' using 1:(log($7)/log(10)) with lines lc 7 notitle, \
     '../data/coexCurveVS_LJ_rc4.dat' using 1:(log($7)/log(10)) with lines lc 8 title 'LJ $r_c=4.0$', \
     '../data/coexCurveLS_LJ_rc4.dat' using 1:(log($7)/log(10)) with lines lc 8 notitle, \
