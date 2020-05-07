set terminal epslatex standalone color size 9cm,7cm
set output 'domega.tex'

#set title 'First derivative'
set title '\Large$kT = 1.2$, \Large$\beta\mu = -3.1$'
set xlabel '\Large$\rho$'
set ylabel '\Large$\beta\partial\omega/\partial\rho$'
set grid linecolor rgb '#B0B0B0'
set key off
set xrange [0:0.6]
set yrange [-0.25:0.25]

plot 'omega_uniform.dat' using 1:3 with lines,\
     'omega_uniform.dat' using 1:(0) with lines linecolor 'black'
