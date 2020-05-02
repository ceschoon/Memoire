set terminal epslatex standalone
set output 'HCP_n.tex'

set title 'HCP Computational Cell (1,1,-2)' offset 0,1
set xlabel '$x$'
set ylabel '$y$'
set zlabel '$z$'

set key off

set xrange[0:1]
set yrange[0:1]
set zrange[0:1]
set xyplane 0

set view 20,135
set view equal xyz
set border 4095

splot 'HCP.dat' using 1:2:3:4 with points pointtype 7 pointsize 10 linecolor variable