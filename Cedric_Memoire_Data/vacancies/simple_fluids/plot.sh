#! /bin/bash

gnuplot   vacancies_simple_fluids.gp
gnuplot   vacancies_simple_fluids_log.gp

pdflatex  vacancies_simple_fluids.tex
pdflatex  vacancies_simple_fluids_log.tex

rm   *.aux   *.log   *.tex  *.eps  *-inc.eps   *-inc-eps-converted-to.pdf



