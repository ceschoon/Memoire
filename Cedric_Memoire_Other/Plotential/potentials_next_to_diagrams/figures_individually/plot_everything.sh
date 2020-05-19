#! /bin/bash

gnuplot   potentials.gp
gnuplot   diagram_WHDF_rc2.gp
gnuplot   diagram_WHDF_rc12.gp

pdflatex  potentials.tex
pdflatex  diagram_WHDF_rc2.tex
pdflatex  diagram_WHDF_rc12.tex

rm   *.aux   *.log   *.tex   *-inc.eps   *-inc-eps-converted-to.pdf

