#! /bin/bash

cd CoexistenceSolid
../../Config.sh
../../dft_make.sh clean
cd ..

cd CoexistenceUniform
../../Config.sh
../../dft_make.sh clean
cd ..

cd DFTAll
../../Config.sh
../../dft_make.sh clean
cd ..

cd DFTFixedkTMuNpoints
../../Config.sh
../../dft_make.sh clean
cd ..

cd DFTMinimizeUniform
../../Config.sh
../../dft_make.sh clean
cd ..

cd Diagram
../../Config.sh
../../dft_make.sh clean
cd ..

cd FindFluidDensitiesFromkTMu
../../Config.sh
../../dft_make.sh clean
cd ..

cd FluidProperties
../../Config.sh
../../dft_make.sh clean
cd ..

cd MinNpoints
../../Config.sh
../../dft_make.sh clean
cd ..

cd QuickTests
./compile.sh
cd ..

cd TestInitialCondition
../../Config.sh
../../dft_make.sh clean
cd ..
