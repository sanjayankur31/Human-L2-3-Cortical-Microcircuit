#!/bin/bash
# Script to clean out generated files
set -x
rm -rf x86_64 arm64
mv *.txt *hoc *.mod *dat *.c *.so *nrn.py *netpyne.py *_eden.py /tmp
