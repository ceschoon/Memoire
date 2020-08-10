#----------------------------------------
set terminal epslatex standalone color size 9cm,7cm
set output 'coexCurveLS.tex'

set title 'Coexistence Liquid-Solid'
set xlabel '$\rho$'
set ylabel '$kT$'
set style data lines
set grid

set key bottom center

plot 'coexCurveLS.dat' using 2:1 title 'liquid', \
     'coexCurveLS.dat' using 3:1 title 'solid'    

