#! /bin/bash

mkdir -p figures

cp TEST_vapour_only/data/*.pdf figures/
cp TEST_vapour_metaliquid/data/*.pdf figures/
cp TEST_coexistence/data/*.pdf figures/
cp TEST_liquid_metavapour/data/*.pdf figures/
cp TEST_liquid_only/data/*.pdf figures/
cp TEST_subcritical/data/*.pdf figures/
cp TEST_critical/data/*.pdf figures/
cp TEST_supercritical/data/*.pdf figures/
cp TEST_supercritical_any/data/*.pdf figures/
