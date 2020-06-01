#! /bin/bash

### This script must be run before running the script "run_diagram"
### I forgot to remove the Cvac field in the report of minimisation results
### The consequence is that it contains values like 1e-312, which cause
### out_of_range errors. The purpose of this script is to remove the Cvac
### field in every coexistenceFS_* files to avoid such error.

zip -r dataLS_backup.zip dataLS
zip -r dataVS_backup.zip dataVS

cd dataLS

for file in $(ls coexistenceFS_*.dat)
do
	cat $file | head -n 9 >  $file.out
	cat $file | tail -n 4 >> $file.out
	mv $file.out $file
done

cd ..

cd dataVS

for file in $(ls coexistenceFS_*.dat)
do
	cat $file | head -n 9 >  $file.out
	cat $file | tail -n 4 >> $file.out
	mv $file.out $file
done

cd ..

