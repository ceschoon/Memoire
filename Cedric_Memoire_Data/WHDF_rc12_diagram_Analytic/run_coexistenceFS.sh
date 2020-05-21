#! /bin/bash

mkdir -p dataFS1
mkdir -p dataFS2

zip -r dataFS1_backup.zip dataFS1
zip -r dataFS2_backup.zip dataFS2

rm dataFS1/*
rm dataFS2/*

cd dataFS1
time ../../../Cedric_Memoire_Code/coexistenceFS/coexistenceFS \
     ../input_files/input_coexistenceFS1.dat &
mv log.dat log1.dat
cd ..

cd dataFS2
time ../../../Cedric_Memoire_Code/coexistenceFS/coexistenceFS \
     ../input_files/input_coexistenceFS2.dat &
mv log.dat log2.dat
cd ..



