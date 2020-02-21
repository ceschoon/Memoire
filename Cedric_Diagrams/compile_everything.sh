#! /bin/bash

cd CoexistenceSolid
#rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..

cd CoexistenceUniform
#rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..

cd DFTAll
#rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..

cd DFTFixedkTMuNpoints
#rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..

cd DFTMinimizeUniform
#rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..

cd Diagram
#rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..

cd FindFluidDensitiesFromkTMu
#rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..

cd FluidProperties
#rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..

cd MinNpoints
#rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..

cd QuickTests
./compile.sh
cd ..

cd TestInitialCondition
#rm -r build debug
../../Config.sh
../../dft_make.sh clean
cd ..
