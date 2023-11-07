#!/usr/bin/env python3
"""
Enter one line description here.

File:

Copyright 2023 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""

from pyneuroml.pynml import generate_nmlgraph

generate_nmlgraph(nml2_file_name="HL23Net_0.5.net.nml", level=3, engine="dot",
                  include_ext_inputs=False)
