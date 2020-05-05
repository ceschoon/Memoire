#! /bin/bash

mkdir -p dataVS
mkdir -p dataLS

cd dataVS
time ../../../Cedric_Memoire_Code/coexistenceFS/coexistenceFS \
     ../input_files/input_coexistenceVS.dat &
cd ..

cd dataLS
time ../../../Cedric_Memoire_Code/coexistenceFS/coexistenceFS \
     ../input_files/input_coexistenceLS.dat &
cd ..


