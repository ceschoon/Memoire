#! /bin/bash

rm -rv build debug
../../Config.sh
../../dft_make.sh

echo "========================="
echo ""
echo "Execute in TEST directory"
echo ""