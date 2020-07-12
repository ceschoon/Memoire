#! /bin/bash

gnuplot   vacancies_gaussian_log.gp
gnuplot   vacancies_with_other_dft_and_simulation.gp
gnuplot   comparison_FCC_HCP.gp

pdflatex  vacancies_gaussian_log.tex
pdflatex  vacancies_with_other_dft_and_simulation.tex
pdflatex  comparison_FCC_HCP.tex

rm   *.aux   *.log   *.tex  *.eps  *-inc.eps   *-inc-eps-converted-to.pdf



