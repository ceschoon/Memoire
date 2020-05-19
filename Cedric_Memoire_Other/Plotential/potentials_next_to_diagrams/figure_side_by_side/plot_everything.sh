#! /bin/bash

gnuplot   multiplot.gp

pdflatex  multiplot.tex

rm   *.aux   *.log   *.tex   *-inc.eps   *-inc-eps-converted-to.pdf

mv multiplot.pdf diagram_variety.pdf