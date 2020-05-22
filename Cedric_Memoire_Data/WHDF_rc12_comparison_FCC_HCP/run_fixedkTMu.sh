#! /bin/bash

mkdir -p fixedkTMu
cd fixedkTMu

### compute HCP 

../../../Cedric_Memoire_Code/fixedkTMu/fixedkTMu \
	../input_files/input_fixedkTMu_HCP1.dat
mv log.dat log_HCP1.dat

../../../Cedric_Memoire_Code/fixedkTMu/fixedkTMu \
	../input_files/input_fixedkTMu_HCP2.dat
mv log.dat log_HCP2.dat

../../../Cedric_Memoire_Code/fixedkTMu/fixedkTMu \
	../input_files/input_fixedkTMu_HCP3.dat
mv log.dat log_HCP3.dat

### compute FCC 

../../../Cedric_Memoire_Code/fixedkTMu/fixedkTMu \
	../input_files/input_fixedkTMu_FCC1.dat
mv log.dat log_FCC1.dat

../../../Cedric_Memoire_Code/fixedkTMu/fixedkTMu \
	../input_files/input_fixedkTMu_FCC2.dat
mv log.dat log_FCC2.dat

../../../Cedric_Memoire_Code/fixedkTMu/fixedkTMu \
	../input_files/input_fixedkTMu_FCC3.dat
mv log.dat log_FCC3.dat



