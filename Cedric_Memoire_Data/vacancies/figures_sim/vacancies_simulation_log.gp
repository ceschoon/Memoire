
set term epslatex standalone color size 9cm,7cm
set output 'vacancies_simulation_log.tex'

set title "Vacancies in simulation (LJ)"
set xlabel "$\\log_{10}kT$"
set ylabel "$\\log_{10}c_{vac}$"
set style data lines
set grid linecolor rgb '#B0B0B0'
#set key bottom right
set key at 2.4,-4.3

plot '../data_dft/coexCurveLS_LJ_rc4_long.dat' using (log($1)/log(10)):(log($7)/log(10)) with lines title 'this work', \
     '../data_sim/out_purohit2018_coex.dat' using (log($4)/log(10)):(log($6)/log(10)):(log($4-$5)/log(10)):(log($4+$5)/log(10)):(log($6-$7)/log(10)):(log($6+$7)/log(10)) with xyerrorbars pointtype 7 linecolor 2 title 'purohit2018'
