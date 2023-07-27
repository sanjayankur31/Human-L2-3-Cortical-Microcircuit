#!/usr/bin/env python3
"""
Test the ProbAMPANMDA synapse

File: NeuroML2/synapses/testProbAMPANMDA.py

Copyright 2023 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""


import sys
import numpy

import neuroml
from neuroml.loaders import read_neuroml2_file
from neuroml.utils import component_factory
from pyneuroml.lems.LEMSSimulation import LEMSSimulation
from pyneuroml.plot import generate_plot
from pyneuroml.pynml import (reload_saved_data, run_lems_with_jneuroml_neuron,
                             write_neuroml2_file)
from matplotlib import pyplot


def main(simulation_id):
    """main method that creates sim and runs it."""
    # a new document
    newdoc = component_factory(neuroml.NeuroMLDocument, id="test_probAMPANMDA")  # type: neuroml.NeuroMLDocument

    # add a cell
    newdoc.add(neuroml.IncludeType(href="passiveCell.cell.nml"))
    # cell = read_neuroml2_file("passiveCell.cell.nml").cells[0]
    # newdoc.add(cell)

    net = newdoc.add(neuroml.Network, id="TestNet", validate=False)

    # Create a population of defined cells and add it to the model
    pop = net.add(neuroml.Population, id="TestPop", component="passiveCell",
                  type="populationList", validate=False)
    pop.add(neuroml.Instance, id=0, location=neuroml.Location(x=0, y=0, z=0))

    # set up a spike array
    spikearray = newdoc.add(neuroml.SpikeArray, id="spikeArray", validate=False)
    idval = 0
    for t in [100, 130, 160]:
        spikearray.add(neuroml.Spike, id=idval, time=f"{t} ms")
        idval += 1

    stimpop = net.add("Population", id="SpikePop", component=spikearray.id, size=1)

    proj = net.add(neuroml.Projection, id="proj", presynaptic_population="SpikePop",
                   postsynaptic_population="TestPop", synapse="probAMPANMDA")
    proj.add(neuroml.ConnectionWD, id=0, pre_cell_id="../SpikePop[0]",
             post_cell_id="../TestPop/0/0", weight="1e-3", delay="0 ms")

    # validate the current document
    newdoc.validate(recursive=True)

    nml_file = "testProbAMPANMDA.net.nml"
    # No longer valid NeuroML because it includes a new Component that the
    # Schema does not know about
    write_neuroml2_file(newdoc, nml_file, validate=True)


    # Create simulation, and record data
    simulation = LEMSSimulation(
        sim_id=simulation_id, duration=500, dt=0.01, simulation_seed=123
    )
    # Include the new Component
    simulation.include_lems_file("ProbAMPANMDA.synapse.xml",
                                 include_included=True)

    simulation.assign_simulation_target(net.id)
    simulation.include_neuroml2_file(nml_file)


    # record other variables
    simulation.create_output_file("output1", f"{simulation_id}.output.dat")
    simulation.add_column_to_output_file("output1", "v", "TestPop/0/0/v")


    simulation.add_column_to_output_file("output1", "i_AMPA", "TestPop/0/0/synapses:probAMPANMDA:0/i_AMPA")
    simulation.add_column_to_output_file("output1", "i_NMDA", "TestPop/0/0/synapses:probAMPANMDA:0/i_NMDA")
    simulation.add_column_to_output_file("output1", "g_AMPA", "TestPop/0/0/synapses:probAMPANMDA:0/g_AMPA")
    simulation.add_column_to_output_file("output1", "g_NMDA", "TestPop/0/0/synapses:probAMPANMDA:0/g_NMDA")
    simulation.add_column_to_output_file("output1", "A_AMPA", "TestPop/0/0/synapses:probAMPANMDA:0/A_AMPA")
    simulation.add_column_to_output_file("output1", "A_NMDA", "TestPop/0/0/synapses:probAMPANMDA:0/A_NMDA")
    simulation.add_column_to_output_file("output1", "B_AMPA", "TestPop/0/0/synapses:probAMPANMDA:0/B_AMPA")
    simulation.add_column_to_output_file("output1", "B_NMDA", "TestPop/0/0/synapses:probAMPANMDA:0/B_NMDA")

    sim_filename = lems_simulation_file = simulation.save_to_file()
    data = run_lems_with_jneuroml_neuron(sim_filename, max_memory="8G", skip_run=False, nogui=True, compile_mods=True, load_saved_data=True)

    return data

def plots(data):
    """Plot bits"""
    print("Generating plots")
    print(f"Data found: {data.keys()}")

    yvalues=[data['TestPop/0/0/v']]
    generate_plot(xvalues=numpy.array([data['t']] * len(yvalues)) * 1000,
                  yvalues=numpy.array(yvalues) * 1000,
                  title="Membrane potential (mV)",
                  labels=["v"], show_plot_already=False)

    yvalues1=[data['TestPop/0/0/synapses:probAMPANMDA:0/A_AMPA'],
              data['TestPop/0/0/synapses:probAMPANMDA:0/A_NMDA'],
              data['TestPop/0/0/synapses:probAMPANMDA:0/B_AMPA'],
              data['TestPop/0/0/synapses:probAMPANMDA:0/B_NMDA']]
    generate_plot(xvalues=[data['t']] * len(yvalues1),
                  yvalues=yvalues1,
                  title="States",
                  labels=[
                          "A_AMPA", "A_NMDA", "B_AMPA", "B_NMDA"
                          ], show_plot_already=False)

    # conductances, multiple by 10e6 to convert to uS to match NEURON mod file
    yvalues2=[data['TestPop/0/0/synapses:probAMPANMDA:0/g_AMPA'],
              data['TestPop/0/0/synapses:probAMPANMDA:0/g_NMDA']
              ]
    generate_plot(xvalues=[data['t']] * len(yvalues2),
                  yvalues=numpy.array(yvalues2),
                  title="Conductances (S)",
                  labels=["g_AMPA", "g_NMDA"], show_plot_already=False)

    pyplot.show()



if __name__ == "__main__":
    simulation_id = "test_probAMPANMDA"

    try:
        data = reload_saved_data(f"LEMS_{simulation_id}.xml")
        print("Data already exists, plotting directly. Remove data files to re-run simulation")
    except OSError:
        print("Data wasn't found. Re-running simulation")
        data = main(simulation_id)
    except Exception:
        print("Data wasn't found. Re-running simulation")
        data = main(simulation_id)

    plots(data)
