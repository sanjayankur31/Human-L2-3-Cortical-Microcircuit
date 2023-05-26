#!/usr/bin/env python3
"""
Get info from generated NEURON files.

This *must* be run after the NEURON files have been generated from the NeuroML.

File: neuroConstruct/getinfo.py

Copyright 2022 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""


import sys
import yaml

from pyneuroml.neuron import morphinfo, getinfo, load_hoc_or_python_file
from neuron import h


if len(sys.argv) != 2:
    print("Error: only takes one argument")
    sys.exit(-1)

cell = sys.argv[1]
load_hoc_or_python_file(f"HL23{cell}.hoc")
h("celsius = 34")
h("objectvar mycell")
h("strdef reference")
h('reference = "acell"')
h(f'mycell = new HL23{cell}(reference, "HL23{cell}", "A cell")')
with open(f"NeuroML-morphinfo-{cell}.yaml", "w") as f:
    retval = morphinfo()
    print(yaml.dump(retval, sort_keys=True, indent=4), file=f, flush=True)
with open(f"NeuroML-info-{cell}.yaml", "w") as f:
    retval = getinfo(h.allsec())
    print(yaml.dump(retval, sort_keys=True, indent=4), file=f, flush=True)
