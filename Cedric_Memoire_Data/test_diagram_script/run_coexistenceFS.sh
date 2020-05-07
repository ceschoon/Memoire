#! /bin/bash

mkdir -p dataVS
mkdir -p dataLS

zip -r dataVS_backup.zip dataVS
zip -r dataLS_backup.zip dataLS

rm dataVS/*
rm dataLS/*

cd dataVS
time ../../../Cedric_Memoire_Code/coexistenceFS/coexistenceFS \
     ../input_files/input_coexistenceVS.dat &
cd ..

cd dataLS
time ../../../Cedric_Memoire_Code/coexistenceFS/coexistenceFS \
     ../input_files/input_coexistenceLS.dat &
cd ..

