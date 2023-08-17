#!/bin/bash
set -ex

./convert_cells_to_neuroml.sh

python postprocess_cells.py -postprocall

rm -rf x86_64 arm64

omv all -V 

python getinfoneuroml.py VIP
python getinfoneuroml.py PV
python getinfoneuroml.py PYR
python getinfoneuroml.py SST

python create_test_network.py -all -nml
python create_test_network.py -net2 -nml


