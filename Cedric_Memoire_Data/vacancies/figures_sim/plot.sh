#! /bin/bash

gnuplot   vacancies_simulation.gp
gnuplot   vacancies_simulation_log.gp

pdflatex  vacancies_simulation.tex
pdflatex  vacancies_simulation_log.tex

rm   *.aux   *.log   *.tex  *.eps  *-inc.eps   *-inc-eps-converted-to.pdf



