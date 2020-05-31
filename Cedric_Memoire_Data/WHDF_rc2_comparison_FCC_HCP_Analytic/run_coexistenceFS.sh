#! /bin/bash

mkdir -p dataVS1
mkdir -p dataLS1
mkdir -p dataLS2

zip -r dataVS1_backup.zip dataVS1
zip -r dataLS1_backup.zip dataLS1
zip -r dataLS2_backup.zip dataLS2

rm dataVS1/*
rm dataLS1/*
rm dataLS2/*

cd dataVS1
time ../../../Cedric_Memoire_Code/coexistenceFS/coexistenceFS \
     ../input_files/input_coexistenceVS1.dat &
cd ..

cd dataLS1
time ../../../Cedric_Memoire_Code/coexistenceFS/coexistenceFS \
     ../input_files/input_coexistenceLS1.dat &
cd ..

cd dataLS2
time ../../../Cedric_Memoire_Code/coexistenceFS/coexistenceFS \
     ../input_files/input_coexistenceLS2.dat &
cd ..

