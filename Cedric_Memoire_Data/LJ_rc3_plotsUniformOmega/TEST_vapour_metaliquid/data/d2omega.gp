set terminal epslatex standalone color size 9cm,7cm
set output 'd2omega.tex'

#set title 'Second derivative'
set title '\Large$kT = 1.2$, \Large$\beta\mu = -3.02$'
set xlabel '\Large$\rho$'
set ylabel '\Large$\beta\partial^2\omega/\partial\rho^2$'
set grid linecolor rgb '#B0B0B0'
set key off
set xrange [0:0.6]
set yrange [-2:3]

plot 'omega_uniform.dat' using 1:4 with lines,\
     'omega_uniform.dat' using 1:(0) with lines linecolor 'black'
