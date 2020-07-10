#set terminal x11
set terminal epslatex standalone size 8cm,5.3cm
set output 'vacancy.tex'

#set title 'FCC Structure' offset 0,-2

#set xlabel '$x$'
#set ylabel '$y$'
set key off

unset xtics
unset ytics
unset border 

#set xrange[0:3.5]
#set yrange[0:2.5]

# borders of the cubic cell
bordercolour = 8
borderwidth = 3
#set arrow from 0,0,0 to 1,0,0 nohead lc bordercolour lw borderwidth

plot 'lattice.dat' using 1:2:4 with points pointtype variable pointsize 5 linecolor 8



