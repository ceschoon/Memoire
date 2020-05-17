
set term epslatex standalone color size 9cm,7cm
set output 'vacancies_simple_fluids_log.tex'

set title "Vacancies along coexistence curves"
set xlabel "$kT$"
set ylabel "$\\log_{10}c_{vac}$"
set style data lines
set grid linecolor rgb '#B0B0B0'
set key bottom right

plot '../data_dft/coexCurveVS_LJ_rc3.dat' using 1:(log($7)/log(10)) with lines linecolor 1 notitle, \
     '../data_dft/coexCurveLS_LJ_rc3.dat' using 1:(log($7)/log(10)) with lines linecolor 1 title 'LJ $r_c=3$', \
     '../data_dft/coexCurveVS_LJ_rc4.dat' using 1:(log($7)/log(10)) with lines linecolor 8 notitle, \
     '../data_dft/coexCurveLS_LJ_rc4.dat' using 1:(log($7)/log(10)) with lines linecolor 8 title 'LJ $r_c=4$', \
     '../data_dft/coexCurveVS_WHDF_rc2.dat' using 1:(log($7)/log(10)) with lines linecolor 2 notitle, \
     '../data_dft/coexCurveLS_WHDF_rc2.dat' using 1:(log($7)/log(10)) with lines linecolor 2 title 'WHDF $r_c=2$', \

