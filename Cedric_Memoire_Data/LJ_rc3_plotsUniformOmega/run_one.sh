#! /bin/bash

suffix=$1

cd TEST_$suffix
../../../Cedric_Memoire_Code/plotsUniformOmega/plotsUniformOmega input.dat
cd data
./plot_everything.sh $suffix
cd ../..
