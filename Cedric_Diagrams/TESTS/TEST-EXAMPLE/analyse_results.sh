#! /bin/bash

# copy the data from all processes into one folder

mkdir -p Data-kt-mu-Npoints-Solid
for PROCESSDIR in $(ls | grep PROCESS)
do
	cp $PROCESSDIR/Data-kt-mu-Npoints-Solid/* Data-kt-mu-Npoints-Solid/
done

# Minimise over Npoints

cp ../../MinNpoints/minNpoints ./
./minNpoints input.dat
mv log.dat log_minNpoints.dat

# Find coexistence

cp ../../CoexistenceSolid/coexistenceSolid ./
./coexistenceSolid input.dat
mv log.dat log_coexistenceSolid.dat

# Draw Diagram

cp ../../Diagram/diagram ./
./diagram