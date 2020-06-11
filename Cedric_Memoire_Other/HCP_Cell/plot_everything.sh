#! /bin/bash

gnuplot   HCP_cell.gp

pdflatex  HCP_cell.tex

rm   *.aux   *.log   *.tex   *-inc.eps   *-inc-eps-converted-to.pdf

