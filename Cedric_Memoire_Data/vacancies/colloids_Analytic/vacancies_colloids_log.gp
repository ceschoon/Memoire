
set term epslatex standalone color size 9cm,7cm
set output 'vacancies_colloids_log.tex'

set title "Vacancies along coexistence curves"
set xlabel "$kT$"
set ylabel "$\\log_{10}c_{vac}$"
set style data lines
set grid linecolor rgb '#B0B0B0'
set key bottom left

set xrange[0.5:1.5]

plot '../data_dft/coexCurveFS_WHDF_rc12_Analytic.dat' using 1:(log($7)/log(10)) with lines linecolor 1 title 'WHDF $r_c=1.2$', \
