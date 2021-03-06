
set term epslatex standalone color size 9cm,7cm
set output 'vacancies_other_dft_with_sim.tex'

set title "Vacancies in other DFT models (LJ)"
set xlabel "$kT$"
set ylabel "$\\log_{10}(c_{vac}/\\rho_S)$"
set style data lines
set grid linecolor rgb '#B0B0B0'
set key bottom left width -3

set xrange[0.6:2.1]
set yrange[-9:-3]

plot '../data_dft/coexCurveLS_LJ_rc4_long.dat' using 1:(log($7/$3)/log(10)) with lines title 'this work', \
     '../data_other_dft/mcrae1990.dat' using 1:2 with linespoints pointtype 7 linecolor 4 title 'mcrae1990', \
     '../data_other_dft/singh2007.dat' using 1:2 with linespoints pointtype 7 linecolor 7 title 'singh2007', \
     '../data_sim/out_purohit2018_coex.dat' using 4:(log(($6)/sqrt($3))/log(10)):($4-$5):($4+$5):(log(($6-$7)/sqrt($3))/log(10)):(log(($6+$7)/sqrt($3))/log(10)) with xyerrorbars pointtype 7 linecolor 2 title 'purohit2018'
