#!/usr/bin/env python3
"""
Enter one line description here.

File:

Copyright 2022 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""

import os

import pyneuroml
from pyneuroml.neuron import export_to_neuroml2

cells = ["HL23PV", "HL23PYR", "HL23SST", "HL23VIP"]


for acell in cells:
    loader_hoc_file = f"{acell}.hoc"
    loader_hoc_file_txt = f"""
    /*load_file("nrngui.hoc")*/
    load_file("stdrun.hoc")

    //=================== creating cell object ===========================
    load_file("import3d.hoc")
    objref newcell

    strdef morphology_file
    morphology_file = "../../L23Net/morphologies/{acell}.swc"

    load_file("../../L23Net/models/biophys_{acell}.hoc")
    load_file("../../L23Net/models/NeuronTemplate.hoc")
    newcell = new NeuronTemplate(morphology_file)
    """

    with open(loader_hoc_file, 'w') as f:
        print(loader_hoc_file_txt, file=f)

    export_to_neuroml2(loader_hoc_file, f"{acell}.cell.nml",
                       includeBiophysicalProperties=False, validate=False)

    os.remove(loader_hoc_file)
