set terminal epslatex standalone color size 9cm,7cm
set output 'omega.tex'

#set title 'Free energy profile'
set title '\Large$kT = 1.35$, \Large$\beta\mu = -2.69$'
set xlabel '$\rho$'
set ylabel '$\beta\omega$'
set grid linecolor rgb '#B0B0B0'
set key off
set xrange [0:0.6]
set yrange [-0.12:0]

plot 'omega_uniform.dat' using 1:2 with lines
