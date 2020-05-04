#! /bin/bash

cd dftSolidComputation
rm -r build debug
../../Config.sh
../../dft_make.sh
cd ..

cd fixedkTMu
rm -r build debug
../../Config.sh
../../dft_make.sh
cd ..

cd surfacePlot
rm -r build debug
../../Config.sh
../../dft_make.sh
cd ..

cd coexistenceFS
rm -r build debug
../../Config.sh
../../dft_make.sh
cd ..

cd coexCurveFS
rm -r build debug
../../Config.sh
../../dft_make.sh
cd ..

cd coexCurveVL
rm -r build debug
../../Config.sh
../../dft_make.sh
cd ..

cd criticalPoint
rm -r build debug
../../Config.sh
../../dft_make.sh
cd ..

cd polyFit
rm -r build debug
../../Config.sh
../../dft_make.sh
cd ..

cd triplePoint
rm -r build debug
../../Config.sh
../../dft_make.sh
cd ..

cd plotsUniformOmega
rm -r build debug
../../Config.sh
../../dft_make.sh
cd ..



