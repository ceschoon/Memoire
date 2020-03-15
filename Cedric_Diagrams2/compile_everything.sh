#! /bin/bash

cd dftSolidComputation
rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..

cd fixedkTMu
rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..

cd surfacePlot
rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..

cd coexistenceFS
rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..

cd coexCurveFS
rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..

cd coexCurveVL
rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..
