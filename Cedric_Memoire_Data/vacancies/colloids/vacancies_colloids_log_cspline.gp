
set term epslatex standalone color size 9cm,7cm
set output 'vacancies_colloids_log_cspline.tex'

set title "Vacancies along coexistence curves \n WHDF $r_c=1.2$"
set xlabel "$kT$"
set ylabel "$\\log_{10}c_{vac}$"
set style data lines
set grid linecolor rgb '#B0B0B0'
set key bottom left width -2

#set label 'WHDF $r_c=1.2$' at 0.23,-3.43
#set label 'WHDF $r_c=1.2$' at 0.4,-2.55

plot '../data_dft/coexCurveFS_WHDF_rc12.dat' using 1:(log($7)/log(10)) with lines linecolor 1 title 'quadratic', \
     '../data_dft/coexCurveFS_WHDF_rc12_cspline.dat' using 1:(log($7)/log(10)) with lines linecolor 2 title 'c. spline', \

