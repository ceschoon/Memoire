#! /bin/bash

gnuplot   potentials.gp

pdflatex   potentials.tex

rm   potentials.aux   potentials.log   potentials.tex   potentials-inc.eps   potentials-inc-eps-converted-to.pdf

