#----------------------------------------
set terminal epslatex standalone color
set output 'surfacePlot.tex'

set title '$kT = 0.791701 \ \beta\mu = -4.76168 \ N_{grid} = 131$'
set xlabel '$\log_{10} c_{vac}$'
set ylabel '$\log_{10}\alpha$'

set key off
set xrange[-8:-1.26058]
set yrange[0:2.91314]

set datafile missing 'failed'
plot 'surfacePlot_kT=7.9170e-01_mu=-4.7617e+00_Ngrid=131.dat' using 1:2:($3/1.0) with image
