#----------------------------------------
set terminal epslatex standalone color
set output 'surfacePlot.tex'

set title '$kT = 0.9 \ \beta\mu = 4.1 \ N_{grid} = 121$'
set xlabel '$\log_{10} c_{vac}$'
set ylabel '$\log_{10}\alpha$'

set key off
set xrange[-8:-1.26058]
set yrange[0:3.86971]

set datafile missing 'failed'
plot 'surfacePlot_kT=9.0000e-01_mu=4.1000e+00_Ngrid=121.dat' using 1:2:($3/1.0) with image
