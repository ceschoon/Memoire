#! /bin/bash

gnuplot   FCC_cell.gp

pdflatex  FCC_cell.tex

rm   *.aux   *.log   *.tex   *-inc.eps   *-inc-eps-converted-to.pdf

