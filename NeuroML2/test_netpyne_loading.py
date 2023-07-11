#!/usr/bin/env python3
"""
Enter one line description here.

File:

Copyright 2023 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""

from netpyne import sim, specs


def load_netpyne():
    """Load a neuroml model into netpyne """
    nml2_file_name = "HL23Net_0.01.net.nml"

    simConfig = specs.SimConfig()
    netParams = sim.importNeuroML2(nml2_file_name, simConfig, simulate=False,
                                   analyze=True, return_net_params_also=True)

    return netParams, simConfig


if __name__ == "__main__":
    np, sc = load_netpyne()
    print(np, sc)
