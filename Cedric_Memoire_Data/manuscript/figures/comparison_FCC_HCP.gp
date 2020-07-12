#----------------------------------------
set terminal epslatex standalone color size 9cm,6cm
set output 'comparison_FCC_HCP.tex'

#set title "FCC vs HCP summary" font ",20"
set xlabel '$kT$'
#set ylabel '$\beta\omega$'
set ylabel '$\beta\omega_{HCP} - \beta\omega_{FCC}$'

set style data lines
set grid
set key bottom right

set yrange[-0.05:0.02]

plot '../data/FCCvsHCP_WHDF_rc2.dat'  using 3:($12-$9)  with points pointtype 7 linecolor 2 title 'WHDF $r_c=2.0$' , \
     '../data/FCCvsHCP_WHDF_rc12.dat' using 3:($12-$9)  with points pointtype 7 linecolor 4 title 'WHDF $r_c=1.2$' , \
