#!/bin/bash

# Copyright 2022 Ankur Sinha
# Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
# File : NeuroML2/scripts/convert_cells_to_neuroml.sh
#
#
# Script to convert cells to NeuroML.


function convert_morphologies() {
    # we must use a shell script to iterate because Neuron needs to be
    # repeatedly quit after each cell.  If it isn't quite and restarted, it
    # keeps previous cells around and so each subsequent export also includes
    # all previous cells.
    #for cell in "HL23PV" "HL23PYR"
    for cell in "HL23PV" "HL23PYR" "HL23SST" "HL23VIP"
    #for cell in "HL23PV"
    do
        python cellmorph2nml.py ${cell}
        #pynml-plotmorph -plane2d "xy" -nogui -saveToFile "${cell}.xy.png" "${cell}.morph.cell.nml"
    done
}

function postprocess () {
    python postprocess_cells.py
}

function clean() {
    echo "Removing: *.hoc *.mod LEMS* *.dat *.nrn.py x86_64 iv*nml"
    rm *.hoc iv*nml
    rm -f *.mod LEMS* *.dat *nrn.py
    rm -rf x86_64
}

function setup () {
    libnrnmechdir="$(dirname $(find . -maxdepth 2 -name "libnrnmech*" ))"
    echo "Removing ${libnrnmechdir}"
    rm -rf "${libnrnmechdir}"
    echo "Compiling mods"
    nrnivmodl mod
}

# Do everything without any arguments
if [ $# -lt 1 ]
then
    setup
    convert_morphologies
    postprocess
    exit 0
fi

# parse options
while getopts "csm" OPTION
do
    case $OPTION in
        c)
            clean
            exit 0
            ;;
        s)
            setup
            exit 0
            ;;
        m)
            clean
            setup
            convert_morphologies
            exit 0
            ;;
        ?)
            echo "Error: unrecognised option"
            exit 1
            ;;
    esac
done
