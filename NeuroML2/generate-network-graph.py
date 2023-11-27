#!/usr/bin/env python3
"""
Enter one line description here.

File:

Copyright 2023 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""

from pyneuroml.pynml import generate_nmlgraph
import matplotlib as mpl


# From pyneuroml.pynml
def generate_matrix_graph(f, level):
    """Generate a connectivity matrix graph.

    Modified version of the code in pyneuroml.pynml
    """
    from neuromllite.MatrixHandler import MatrixHandler
    from neuroml.hdf5.NeuroMLXMLParser import NeuroMLXMLParser

    handler = MatrixHandler(level=level, nl_network=None)
    currParser = NeuroMLXMLParser(handler)
    currParser.parse(f)
    handler.finalise_document()


if __name__ == "__main__":
    filename = "HL23Net_0.5.net.nml"
    mpl.rcParams['axes.titlesize'] = 30
    mpl.rcParams['axes.labelsize'] = 30
    mpl.rcParams['xtick.labelsize'] = 20
    mpl.rcParams['ytick.labelsize'] = 20
    generate_matrix_graph(filename, 1)

    generate_nmlgraph(nml2_file_name=filename, level=3, engine="dot",
                      include_ext_inputs=False)
