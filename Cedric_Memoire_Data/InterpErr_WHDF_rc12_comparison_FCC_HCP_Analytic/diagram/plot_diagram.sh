#! /bin/bash

gnuplot   diagram.gp

pdflatex  diagram.tex

rm   diagram.aux   diagram.log   diagram.tex   diagram-inc.eps   diagram-inc-eps-converted-to.pdf

mv diagram.pdf diagram_FCC_HCP.pdf
