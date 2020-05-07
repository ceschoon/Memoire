#----------------------------------------
set terminal epslatex standalone color size 9cm,7cm
set output 'coexCurveVL.tex'

set title 'Coexistence Vapour-Liquid'
set xlabel '$\rho$'
set ylabel '$kT$'
set style data lines
set grid

set key bottom left

plot 'coexCurveVL.dat' using 2:1 title 'vapour', \
     'coexCurveVL.dat' using 3:1 title 'liquid'    

