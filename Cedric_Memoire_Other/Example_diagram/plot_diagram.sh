#! /bin/bash

gnuplot   diagram_example.gp

pdflatex  diagram_example.tex

rm   diagram*.aux   diagram*.log   diagram*.tex   diagram*-inc.eps   diagram*-inc-eps-converted-to.pdf


