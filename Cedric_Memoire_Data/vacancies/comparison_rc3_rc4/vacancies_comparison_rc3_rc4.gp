
set term epslatex standalone color size 9cm,7cm
set output 'vacancies_comparison_rc3_rc4.tex'

set title "Effect of the cut-off (LJ)"
set xlabel "$kT$"
set ylabel "$c_{vac} \\times 10^{-4}$"
set style data lines
set grid linecolor rgb '#B0B0B0'
set key top right

#set xrange[0.6:2.1]
set xrange[0.6:1.6]

plot '../data_dft/coexCurveLS_LJ_rc4.dat' using 1:($7*1e4) with lines title '$r_c=4$', \
     '../data_dft/coexCurveLS_LJ_rc3.dat' using 1:($7*1e4) with lines title '$r_c=3$'