#----------------------------------------
set terminal epslatex standalone color size 9cm,7cm
set output 'plot.tex'

set title "FCC vs HCP (WHDF $r_c=1.2$)" font ",20"
set xlabel '$kT$'
#set ylabel '$\beta\omega$'
set ylabel '$\beta\omega_{HCP} - \beta\omega_{FCC}$'

set style data lines
set grid
set key bottom left

#plot 'data.dat' using 3:9  with points pointtype 7 linecolor 1 title 'FCC', \
#     'data.dat' using 3:12 with points pointtype 7 linecolor 2 title 'HCP', \

plot 'data.dat' using 3:($12-$9)  with points pointtype 7 linecolor 1 notitle

