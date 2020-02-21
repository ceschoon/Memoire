#! /bin/bash

g++ -std=c++11 -o readColumnVectorFromFile readColumnVectorFromFile.cpp ../utilities.cpp -lgsl -lgslcblas -lm
g++ -std=c++11 -o test_zeroFromData        test_zeroFromData.cpp        ../utilities.cpp -lgsl -lgslcblas -lm
g++ -std=c++11 -o test_minFromDataParabola test_minFromDataParabola.cpp ../utilities.cpp -lgsl -lgslcblas -lm
