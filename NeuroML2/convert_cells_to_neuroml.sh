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
    #for cell in "HL23PV" "HL23PYR" "HL23SST" "HL23VIP"
    for cell in "HL23PV"
    do
        python cellmorph2nml.py ${cell}
        #pynml-plotmorph -plane2d "xy" -nogui -saveToFile "${cell}.xy.png" "${cell}.morph.cell.nml"
    done
}

function addbiophy () {
    python addcellbiophysics.py
}

convert_morphologies
addbiophy
