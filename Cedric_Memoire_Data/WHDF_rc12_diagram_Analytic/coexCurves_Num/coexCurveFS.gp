#----------------------------------------
set terminal epslatex standalone color size 9cm,7cm
set output 'coexCurveFS.tex'

set title 'Coexistence Fluid-Solid'
set xlabel '$\rho$'
set ylabel '$kT$'
set style data lines
set grid

set key bottom center

plot 'coexCurveFS.dat' using 2:1 title 'fluid', \
     'coexCurveFS.dat' using 3:1 title 'solid'    

