set terminal epslatex standalone color size 6cm,5cm
set output 'potentials.tex'

#set title 'Potential shapes'
set xlabel '\Large$r$'
set ylabel '\Large$v(r)$'
#set grid linecolor rgb '#B0B0B0'
set key off
#set xrange [-0.2:2.2]
set xrange [0.8:2.2]
set yrange [-1.2:1.2]
unset xtics
unset ytics
set xzeroaxis linetype 1 linecolor 8
set yzeroaxis linetype 1 linecolor 8

plot "../data/potentials.dat" using 1:3 with lines linecolor 1 notitle ,\
     "../data/potentials.dat" using 1:4 with lines linecolor 2 notitle
