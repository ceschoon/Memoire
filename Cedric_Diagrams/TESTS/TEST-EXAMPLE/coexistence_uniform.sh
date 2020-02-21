#! /bin/bash

# Tip: move weights to coexistence directory to avoid regenereting them

cp ../../CoexistenceUniform/coexistenceUniform ./
cd CoexistenceUniform
time nice ../coexistenceUniform input.dat > /dev/null &