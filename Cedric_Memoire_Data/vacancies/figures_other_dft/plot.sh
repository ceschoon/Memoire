#! /bin/bash

gnuplot   vacancies_other_dft.gp
gnuplot   vacancies_other_dft_with_sim.gp

pdflatex  vacancies_other_dft.tex
pdflatex  vacancies_other_dft_with_sim.tex

rm   *.aux   *.log   *.tex  *.eps  *-inc.eps   *-inc-eps-converted-to.pdf



