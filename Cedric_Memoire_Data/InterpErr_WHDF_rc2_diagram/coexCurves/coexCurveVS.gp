#----------------------------------------
set terminal epslatex standalone color size 9cm,7cm
set output 'coexCurveVS.tex'

set title 'Coexistence Vapour-Solid'
set xlabel '$\rho$'
set ylabel '$kT$'
set style data lines
set grid

set key bottom center

plot 'coexCurveVS.dat' using 2:1 title 'vapour', \
     'coexCurveVS.dat' using 3:1 title 'solid'    

