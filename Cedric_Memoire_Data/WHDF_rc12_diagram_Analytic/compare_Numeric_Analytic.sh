#! /bin/bash

### fetch data from numeric weights

cp -r ../WHDF_rc12_diagram/coexCurves coexCurves_Num

### draw phase diagram

cd diagram_comparison
./plot_diagram.sh
cd ..

