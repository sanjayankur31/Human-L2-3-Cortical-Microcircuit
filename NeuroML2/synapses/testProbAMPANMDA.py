#!/usr/bin/env python3
"""
Test the ProbAMPANMDA synapse

File: NeuroML2/synapses/testProbAMPANMDA.py

Copyright 2023 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""


import sys
import neuroml
from neuroml.utils import component_factory
from neuroml.loaders import read_neuroml2_file
from pyneuroml.pynml import write_neuroml2_file
from pyneuroml.lems.LEMSSimulation import LEMSSimulation

# a new document
newdoc = component_factory(neuroml.NeuroMLDocument, id="test_probAMPANMDA")

# include the new synapse type
newdoc.add(neuroml.IncludeType, href="ProbAMPANMDA.synapsedef.nml")
newdoc.add(neuroml.IncludeType, href="ProbAMPANMDA.synapse.nml")

# the syn component that I have
# syncomp = read_neuroml2_file("ProbAMPANMDA.synapse.nml")
# does not recognise the new component type, so it won't show the created
# component
# print(syncomp)

# add a cell
izh0 = newdoc.add(neuroml.Izhikevich2007Cell,
                  id="izh2007RS0", v0="-60mV", C="100pF", k="0.7nS_per_mV",
                  vr="-60mV", vt="-40mV", vpeak="35mV", a="0.03per_ms", b="-2nS",
                  c="-50.0mV", d="100pA"
                  )

net = newdoc.add(neuroml.Network, id="TestNet", validate=False)

# Create a population of defined cells and add it to the model
pop = net.add(neuroml.Population, id="TestPop", component=izh0.id, size=1)

# set up a spike array
spikearray = newdoc.add(neuroml.SpikeArray, id="spikeArray", validate=False)
idval = 0
for t in range(30, 500, 10):
    spikearray.add(neuroml.Spike, id=idval, time=f"{t} ms")
    idval += 1

stimpop = net.add("Population", id="SpikePop", component=spikearray.id, size=1)

proj = net.add(neuroml.Projection, id="proj", presynaptic_population="SpikePop",
               postsynaptic_population="TestPop", synapse="probAMPANMDASyn")
proj.add(neuroml.Connection, id=0, pre_cell_id="../SpikePop[0]",
         post_cell_id="../TestPop[0]")

nml_file = "testProbAMPANMDA.net.nml"
write_neuroml2_file(newdoc, nml_file, validate=False)


# Create simulation, and record data
simulation_id = "test_probAMPANMDA"
simulation = LEMSSimulation(
    sim_id=simulation_id, duration=1000, dt=0.1, simulation_seed=123
)
simulation.assign_simulation_target(net.id)
simulation.include_neuroml2_file(nml_file)

# record spikes
simulation.create_event_output_file("output0", f"{simulation_id}.spikes.dat")
simulation.add_selection_to_event_output_file("output0", 0, "TestPop[0]", "spike")

# record other variables
simulation.create_output_file("output1", f"{simulation_id}.output.dat")
simulation.add_column_to_output_file("output1", "v", "TestPop[0]/v")

lems_simulation_file = simulation.save_to_file()
