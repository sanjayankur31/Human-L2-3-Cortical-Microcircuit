#!/usr/bin/env python3
"""
Get info from NEURON model

File: NEURON/getinfo.py 

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
load_hoc_or_python_file(f"test_HL23{cell}.hoc")

with open(f"NEURON-morphinfo-{cell}.yaml", "w") as f:
    retval = morphinfo()
    print(yaml.dump(retval, sort_keys=True, indent=4), file=f, flush=True)
with open(f"NEURON-info-{cell}.yaml", "w") as f:
    retval = getinfo(h.allsec())
    print(yaml.dump(retval, sort_keys=True, indent=4), file=f, flush=True)
