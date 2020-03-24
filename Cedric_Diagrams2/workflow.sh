############################################################################
##                                                                        ##
##           This is a script that shows how to use this code             ##
##           for tracing diagrams. It is not meant to be run              ##
##           at once but rather to show what does what and how            ##
##           to navigate in the folders. The code is not the              ##
##           most stable so it is better to run each command              ##
##           seperataly in a terminal to see if everything goes           ##
##           the right way.                                               ##
##                                                                        ##
############################################################################

##########
## Comment this line to run the entire script (not recommended)
exit

##########
## Make sure that everything is compiled

./compile_everything.sh

########## 
## coexCurveVL: generate data and curve for VL coexistence

cd coexCurveVL/TEST/
../coexCurveVL input.dat
cd ../..

##########
## polyFit: fit the raw data using polynomials

cp coexCurveVL/TEST/coexCurveVL.dat polyFit/TEST/data/
cd polyFit/TEST
../polyFit input_VL_rhoV.dat
../polyFit input_VL_rhoL.dat
../polyFit input_VL_mu.dat
../polyFit input_VL_omega.dat
cd ../..

##########
## criticalPoint: compute properties at the critical point

##########
## coexistenceVS: generate data for VS coexistence

cd coexistenceVS/TEST/
../coexistenceVS input.dat
cd ../..

##########
## coexistenceLS: generate data for LS coexistence

cd coexistenceLS/TEST/
../coexistenceLS input.dat
cd ../..

##########
## coexistenceFS: generate data for FS coexistence

cd coexistenceFS/TEST/
../coexistenceFS input.dat
cd ../..

##########
## coexCurveVS: generate curve for VS coexistence (from existing data)

rm coexCurveVS/TEST/data/*
cp coexistenceVS/TEST/data/* coexCurveVS/TEST/data/
cd coexCurveVS/TEST
../coexCurveVS input.dat
cd ../..

##########
## coexCurveLS: generate curve for LS coexistence (from existing data)

rm coexCurveLS/TEST/data/*
cp coexistenceLS/TEST/data/* coexCurveLS/TEST/data/
cd coexCurveLS/TEST
../coexCurveLS input.dat
cd ../..

##########
## coexCurveFS: generate curve for FS coexistence (from existing data)

rm coexCurveFS/TEST/data/*
cp coexistenceFS/TEST/data/* coexCurveFS/TEST/data/
cd coexCurveFS/TEST
../coexCurveFS input.dat
cd ../..

##########
## polyFit: fit the raw data using polynomials

cp coexCurveVS/TEST/coexCurveVS.dat polyFit/TEST/data/
cd polyFit/TEST
../polyFit input_VS_rhoV.dat
../polyFit input_VS_rhoS.dat
../polyFit input_VS_mu.dat
../polyFit input_VS_omega.dat
cd ../..

##########
## polyFit: fit the raw data using polynomials

cp coexCurveLS/TEST/coexCurveLS.dat polyFit/TEST/data/
cd polyFit/TEST
../polyFit input_LS_rhoL.dat
../polyFit input_LS_rhoS.dat
../polyFit input_LS_mu.dat
../polyFit input_LS_omega.dat
cd ../..

##########
## polyFit: fit the raw data using polynomials

cp coexCurveFS/TEST/coexCurveFS.dat polyFit/TEST/data/
cd polyFit/TEST
../polyFit input_FS_rhoF.dat
../polyFit input_FS_rhoS.dat
../polyFit input_FS_mu.dat
../polyFit input_FS_omega.dat
cd ../..

##########
## triplePoint: compute properties at the triple point

cp polyFit/TEST/data/* triplePoint/TEST/data/*
cd triplePoint/TEST/
../triplePoint input.dat
cd ../..

##########
## drawDiagram: generate a plot of the full phase diagram

cp coexCurveVL/TEST/coexCurveVL.dat drawDiagram/
cp coexCurveVS/TEST/coexCurveVS.dat drawDiagram/
cp coexCurveLS/TEST/coexCurveLS.dat drawDiagram/
cp coexCurveFS/TEST/coexCurveFS.dat drawDiagram/
cd drawdiagram/
gnuplot plot_diagram
cd ..

