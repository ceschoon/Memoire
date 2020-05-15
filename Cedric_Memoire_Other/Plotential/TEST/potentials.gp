set terminal epslatex standalone color size 12cm,9cm
set output 'potentials.tex'

#set title 'Potential shapes'
set xlabel '\Large$r/\sigma$'
set ylabel '\Large$v(r)/\varepsilon$'
set grid linecolor rgb '#B0B0B0'
set key top right
set xrange [0.9:2.1]
set yrange [-1.2:1.2]

plot "potentials.dat" using 1:2 with lines linecolor 8title "LJ $r_c=3$" ,\
     "potentials.dat" using 1:3 with lines linecolor 1title "WHDF $r_c=2$" ,\
     "potentials.dat" using 1:4 with lines linecolor 2title "WHDF $r_c=1.2$" ,\
