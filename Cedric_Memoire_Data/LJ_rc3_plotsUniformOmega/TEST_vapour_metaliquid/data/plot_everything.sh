#! /bin/bash

gnuplot   omega.gp
gnuplot  domega.gp
gnuplot d2omega.gp

pdflatex   omega.tex
pdflatex  domega.tex
pdflatex d2omega.tex

rm   omega.aux   omega.log   omega.tex   omega-inc.eps   omega-inc-eps-converted-to.pdf
rm  domega.aux  domega.log  domega.tex  domega-inc.eps  domega-inc-eps-converted-to.pdf
rm d2omega.aux d2omega.log d2omega.tex d2omega-inc.eps d2omega-inc-eps-converted-to.pdf

if [ $# -gt 0 ]
then
mv   omega.pdf   omega_$1.pdf
mv  domega.pdf  domega_$1.pdf
mv d2omega.pdf d2omega_$1.pdf
fi

