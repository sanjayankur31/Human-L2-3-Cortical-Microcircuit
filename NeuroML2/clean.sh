#!/bin/bash
# Script to clean out generated files
set -ex
rm -rf x86_64 arm64
atmpdir=$(mktemp -d -t nml_)
mv *hoc *.mod *dat *.c *.so *nrn.py *netpyne.py *_eden.py Sim_IC* ${atmpdir}
