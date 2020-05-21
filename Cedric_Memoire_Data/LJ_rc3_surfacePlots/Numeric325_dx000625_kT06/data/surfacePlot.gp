#----------------------------------------
set terminal epslatex standalone color
set output 'surfacePlot.tex'

set title '$kT = 0.6 \ \beta\mu = -7.6 \ N_{grid} = 256$'
set xlabel '$\log_{10} c_{vac}$'
set ylabel '$\log_{10}\alpha$'

set key off
set xrange[-8:-1.26058]
set yrange[0:2.91314]

set datafile missing 'failed'
plot 'surfacePlot_kT=6.0000e-01_mu=-7.6000e+00_Ngrid=256.dat' using 1:2:($3/1.0) with image
