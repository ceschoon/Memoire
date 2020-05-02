#! /bin/bash

gnuplot   FCC.gp
gnuplot   FCC_111.gp
gnuplot   FCC_111_wide.gp

gnuplot   HCP.gp
gnuplot   HCP_n.gp
gnuplot   HCP_n_wide.gp

pdflatex  FCC.tex
pdflatex  FCC_111.tex
pdflatex  FCC_111_wide.tex

pdflatex  HCP.tex
pdflatex  HCP_n.tex
pdflatex  HCP_n_wide.tex

rm   *.aux   *.log   *.tex   *-inc.eps   *-inc-eps-converted-to.pdf

