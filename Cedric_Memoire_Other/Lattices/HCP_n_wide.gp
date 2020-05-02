set terminal epslatex standalone
set output 'HCP_n_wide.tex'

set title 'HCP Lattice (1,-1,0)'
set xlabel '$x$'
set ylabel '$y$'
set zlabel '$z$'

set key off

#set xrange[0:1]
#set yrange[0:1]
#set zrange[0:1]
set xyplane 0

set view 20,135
set view equal xyz

splot 'HCP.dat' using 1:2:3:4 with points pointtype 7 pointsize 2 linecolor variable
