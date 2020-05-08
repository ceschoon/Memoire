#! /bin/bash

gnuplot   diagram_LJ_rc3.gp
gnuplot   diagram_WHDF_rc2.gp
gnuplot   diagram_WHDF_rc12.gp

pdflatex  diagram_LJ_rc3.tex
pdflatex  diagram_WHDF_rc2.tex
pdflatex  diagram_WHDF_rc12.tex

rm   diagram*.aux   diagram*.log   diagram*.tex   diagram*-inc.eps   diagram*-inc-eps-converted-to.pdf

mv diagram_LJ_rc3.pdf    diagram_LJ_rc3_with_sim.pdf
mv diagram_WHDF_rc2.pdf  diagram_WHDF_rc2_with_sim.pdf
mv diagram_WHDF_rc12.pdf diagram_WHDF_rc12_with_sim.pdf

