#----------------------------------------
set terminal epslatex standalone color
set output 'surfacePlot.tex'

set title '$kT = 1.5 \ \beta\mu = 4.6 \ N_{grid} = 128$'
set xlabel '$\log_{10} c_{vac}$'
set ylabel '$\log_{10}\alpha$'

set key off
set xrange[-8:-1.26058]
set yrange[0:2.91314]

set datafile missing 'failed'
plot 'surfacePlot_kT=1.5000e+00_mu=4.6000e+00_Ngrid=128.dat' using 1:2:($3/1.0) with image
