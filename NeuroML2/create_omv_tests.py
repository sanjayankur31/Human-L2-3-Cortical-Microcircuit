#!/usr/bin/env python3
"""
Create OMV test templates

File: create_omv_tests.py

Copyright 2023 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""


import textwrap


# cellnames = ["HL23PV", "HL23PYR", "HL23SST", "HL23VIP"]
cellnames = ["HL23PYR", "HL23SST", "HL23VIP"]
sim_engines = {
    "jnmleden": "jNeuroML_EDEN",
    "jnmlnetpyne": "jNeuroML_NetPyNE",
    "jnmlneuron": "jNeuroML_NEURON",
}

for cell in cellnames:
    cell_suffix = cell[4:].lower()

    for fn, eng in sim_engines.items():
        omv_filename = f".test.{cell_suffix}.{fn}.omt"
        print(f"Creating {omv_filename}")
        with open(omv_filename, 'w') as omv_file:
            print(textwrap.dedent(f"""
# Script for running automated tests on OSB, see https://github.com/OpenSourceBrain/osb-model-validation

target: LEMS_{cell}_sim.xml
engine: {eng}
mep: tests/.test.{cell_suffix}.mep
experiments:
  step{cell[4:]}:
    observables:
      spike times:
        file:
          path: {cell}_net.dat
          columns: [0,1]
          scaling: [1000,1000]
        spike detection:
          method: threshold
            """),
                  file=omv_file)
        # tolerance: 0.0
