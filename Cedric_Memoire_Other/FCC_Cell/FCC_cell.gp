#set terminal x11
set terminal epslatex standalone size 8cm,8cm
set output 'FCC_cell.tex'

#set title 'FCC Structure' offset 0,-2

#set xlabel '$x$'
#set ylabel '$y$'
#set zlabel '$z$'

unset xtics
unset ytics
unset ztics

set key off

#set xrange[-0.5:1.5]
#set yrange[-0.5:2]
#set zrange[0:2]
set xyplane 0

set view 70,15
set view equal xyz
unset border 
#set border 4095



# borders of the cubic cell
bordercolour = 2
borderwidth = 3

# connect lower square
set arrow from 0,0,0 to 1,0,0 nohead lc bordercolour lw borderwidth
set arrow from 1,0,0 to 1,1,0 nohead lc bordercolour lw borderwidth
set arrow from 1,1,0 to 0,1,0 nohead lc bordercolour lw borderwidth
set arrow from 0,1,0 to 0,0,0 nohead lc bordercolour lw borderwidth

# connect upper square
set arrow from 0,0,1 to 1,0,1 nohead lc bordercolour lw borderwidth
set arrow from 1,0,1 to 1,1,1 nohead lc bordercolour lw borderwidth
set arrow from 1,1,1 to 0,1,1 nohead lc bordercolour lw borderwidth
set arrow from 0,1,1 to 0,0,1 nohead lc bordercolour lw borderwidth

# connect the two squares
set arrow from 0,0,0 to 0,0,1 nohead lc bordercolour lw borderwidth
set arrow from 1,0,0 to 1,0,1 nohead lc bordercolour lw borderwidth
set arrow from 1,1,0 to 1,1,1 nohead lc bordercolour lw borderwidth
set arrow from 0,1,0 to 0,1,1 nohead lc bordercolour lw borderwidth




splot 'FCC_cell.dat' using 1:2:3:4 with points pointtype 7 pointsize 3 linecolor variable



