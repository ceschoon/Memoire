#set terminal x11
set terminal epslatex standalone size 8cm,8cm
set output 'HCP_cell.tex'

#set title 'HCP Structure' offset 0,-2

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

set view 75,10
set view equal xyz
unset border 
#set border 4095



# connect lower hexagon
set arrow from 1,0,0 to 1.5,0.866,0 nohead
set arrow from 1.5,0.866,0 to 1,1.732,0 nohead
set arrow from 1,1.732,0 to 0,1.732,0 nohead
set arrow from 0,1.732,0 to -0.5,0.866,0 nohead
set arrow from -0.5,0.866,0 to 0,0,0 nohead
set arrow from 0,0,0 to 1,0,0 nohead

# connect upper hexagon
set arrow from 1,0,1.633 to 1.5,0.866,1.633 nohead
set arrow from 1.5,0.866,1.633 to 1,1.732,1.633 nohead
set arrow from 1,1.732,1.633 to 0,1.732,1.633 nohead
set arrow from 0,1.732,1.633 to -0.5,0.866,1.633 nohead
set arrow from -0.5,0.866,1.633 to 0,0,1.633 nohead
set arrow from 0,0,1.633 to 1,0,1.633 nohead

# connect the two hexagons
set arrow from 1,0,0 to 1,0,1.633 nohead
set arrow from 1.5,0.866,0 to 1.5,0.866,1.633 nohead
set arrow from 1,1.732,0 to 1,1.732,1.633 nohead
set arrow from 0,1.732,0 to 0,1.732,1.633 nohead
set arrow from -0.5,0.866,0 to -0.5,0.866,1.633 nohead
set arrow from 0,0,0 to 0,0,1.633nohead



# borders of the rectangular cell

bordercolour = 2
borderwidth = 3

set arrow from 0,0,0 to 1,0,0 nohead lc bordercolour lw borderwidth
set arrow from 1,0,0 to 1,1.732,0 nohead lc bordercolour lw borderwidth
set arrow from 1,1.732,0 to 0,1.732,0 nohead lc bordercolour lw borderwidth
set arrow from 0,1.732,0 to 0,0,0 nohead lc bordercolour lw borderwidth

set arrow from 0,0,1.633 to 1,0,1.633 nohead lc bordercolour lw borderwidth
set arrow from 1,0,1.633 to 1,1.732,1.633 nohead lc bordercolour lw borderwidth
set arrow from 1,1.732,1.633 to 0,1.732,1.633 nohead lc bordercolour lw borderwidth
set arrow from 0,1.732,1.633 to 0,0,1.633 nohead lc bordercolour lw borderwidth

set arrow from 0,0,0 to 0,0,1.633 nohead lc bordercolour lw borderwidth
set arrow from 1,0,0 to 1,0,1.633 nohead lc bordercolour lw borderwidth
set arrow from 1,1.732,0 to 1,1.732,1.633 nohead lc bordercolour lw borderwidth
set arrow from 0,1.732,0 to 0,1.732,1.633 nohead lc bordercolour lw borderwidth



splot 'HCP_cell.dat' using 1:2:3:4 with points pointtype 7 pointsize 2 linecolor variable



