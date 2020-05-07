#! /bin/bash

for dir in $(./run_list.sh)
do
	cd $dir
	nice -n 10 time ../../../Cedric_Memoire_Code/surfacePlot/surfacePlot input.dat &
	cd ..
done
