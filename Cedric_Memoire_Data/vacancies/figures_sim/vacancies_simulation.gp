
set term epslatex standalone color size 9cm,7cm
set output 'vacancies_simulation.tex'

set title "Vacancies in simulation (LJ)"
set xlabel "$kT$"
set ylabel "$c_{vac} \\times 10^{-4}$"
set style data lines
set grid linecolor rgb '#B0B0B0'
#set key at 3.7,1.6 
set key top right

set xrange[0.6:3.8]
set yrange[0:4]

plot '../data_dft/coexCurveLS_LJ_rc4_long.dat' using 1:($7*1e4) with lines title 'this work', \
     '../data_sim/out_purohit2018_coex.dat' using 4:($6*1e4):($4-$5):($4+$5):($6*1e4-$7*1e4):($6*1e4+$7*1e4) with xyerrorbars pointtype 7 linecolor 2 title 'purohit2018'