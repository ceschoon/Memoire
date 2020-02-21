#! /bin/bash

# Distributed over multiple processes to go faster

cp ../../DFTAll/dftAll_gaussian ./

for PROCESSDIR in $(ls | grep PROCESS)
do
	cd $PROCESSDIR
	time nice ../dftAll_gaussian input.dat > /dev/null &
	cd ..
done
