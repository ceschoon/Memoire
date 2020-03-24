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

cd coexistenceVS
rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..

cd coexistenceLS
rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..

cd coexCurveFS
rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..

cd coexCurveVS
rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..

cd coexCurveLS
rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..

cd coexCurveVL
rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..

cd polyFit
rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..

cd triplePoint
rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..


