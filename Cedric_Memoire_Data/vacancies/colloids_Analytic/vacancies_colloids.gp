
set term epslatex standalone color size 9cm,7cm
set output 'vacancies_colloids.tex'

set title "Vacancies along coexistence curves"
set xlabel "$kT$"
set ylabel "$c_{vac} \\times 10^{-4}$"
set style data lines
set grid linecolor rgb '#B0B0B0'
#set key at 3.7,1.6 
set key top right

set xrange[0.5:1.5]
set yrange[0:3]

plot '../data_dft/coexCurveFS_WHDF_rc12_Analytic.dat' using 1:($7*1e4) with lines linecolor 1 title 'WHDF $r_c=1.2$', \