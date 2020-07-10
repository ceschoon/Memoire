#! /bin/bash

g++ generate_lattice.cpp
./a.out
rm a.out

gnuplot   vacancy.gp

pdflatex  vacancy.tex

rm   *.aux   *.log   *.tex   *-inc.eps   *-inc-eps-converted-to.pdf

convert -density 1000 vacancy.pdf vacancy.png
#convert vacancy.png -crop 1325x875+950+400 vacancy_cropped.png
convert vacancy.png -crop 1360x950+920+430 vacancy_cropped.png

