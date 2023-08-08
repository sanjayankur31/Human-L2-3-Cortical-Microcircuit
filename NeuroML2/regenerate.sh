#!/bin/bash
set -ex

python create_test_network.py -all -nml
python create_test_network.py -net2 -nml


