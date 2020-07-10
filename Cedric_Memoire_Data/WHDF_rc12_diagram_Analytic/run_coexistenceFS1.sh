#! /bin/bash

mkdir -p dataFS1

zip -r dataFS1_backup.zip dataFS1

rm dataFS1/*

cd dataFS1
time ../../../Cedric_Memoire_Code/coexistenceFS/coexistenceFS \
     ../input_files/input_coexistenceFS1.dat &
mv log.dat log1.dat
cd ..


