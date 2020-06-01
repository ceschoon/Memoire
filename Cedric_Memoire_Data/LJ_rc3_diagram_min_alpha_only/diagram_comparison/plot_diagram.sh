#! /bin/bash

gnuplot   diagram.gp

pdflatex  diagram.tex

rm   *.aux   *.log   *.tex   *-inc.eps   *-inc-eps-converted-to.pdf

mv diagram.pdf diagram_LJ_rc3_no_min_cvac.pdf

