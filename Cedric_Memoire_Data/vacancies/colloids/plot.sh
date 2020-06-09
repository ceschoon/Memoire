#! /bin/bash

gnuplot   vacancies_colloids.gp
gnuplot   vacancies_colloids_log.gp
gnuplot   vacancies_colloids_log_cspline.gp

pdflatex  vacancies_colloids.tex
pdflatex  vacancies_colloids_log.tex
pdflatex  vacancies_colloids_log_cspline.tex

rm   *.aux   *.log   *.tex  *.eps  *-inc.eps   *-inc-eps-converted-to.pdf



