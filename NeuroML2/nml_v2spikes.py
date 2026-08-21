#!/usr/bin/env python3
"""
Get and plot spikes from membrane potential data, using neo and elephant.

File: NeuroML2/nml_v2spikes.py

Copyright 2025 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""

import sys

import elephant
import quantities
from lems.model.model import Model
from neo.io import AsciiSignalIO
from pyneuroml.io import read_neuroml2_file
from pyneuroml.plot.PlotSpikes import plot_spikes
from pyneuroml.utils.units import get_value_in_si

if len(sys.argv) < 2:
    print("Script takes one argument: name of the LEMS simulation file")
    sys.exit(1)


print(f"Processing simulation file: {sys.argv[1]}")

model = Model(include_includes=True, fail_on_missing_includes=False)
model.import_from_file(sys.argv[1])

# Note: could have also just parsed the Components here too
network_file = ""
for inc in model.included_files:
    if inc.endswith(".net.nml"):
        network_file = inc

network_doc = read_neuroml2_file(network_file)
populations = network_doc.networks[0].populations

cells = {}
total_cells = 0
for pop in populations:
    cell_type = pop.id.replace("HL23", "").replace("_pop", "")
    pop_size = len(pop.instances)
    cells[cell_type] = pop_size
    total_cells += pop_size


sampling_rate = 1 * quantities.Hz
data_filename = ""

simulation_components = model.get_component_list("Sim")
for com_id, comp in simulation_components.items():
    if comp.type == "Simulation":
        params = comp.parameters
        sampling_rate = (1 / get_value_in_si(params["step"])) * quantities.Hz

        for ch in comp.children:
            if ch.type == "OutputFile":
                params = ch.parameters
                data_filename = params["fileName"]

print(f"Cells in model are: {cells}")
print(f"Samping rate is: {sampling_rate}")
print(f"Data file is: {data_filename}")

assert len(network_file) > 0
assert len(data_filename) > 0

# 200 colums, but skip the 0th time column, and the spurious last column
cols = list(range(1, total_cells + 1))

# sampling rate is 1/dt
memb_pot_data = AsciiSignalIO(
    data_filename,
    timecolumn=None,
    sampling_rate=sampling_rate,
    units=quantities.volt,
    time_units=quantities.second,
    t_start=0 * quantities.second,
    signal_group_mode="all-in-one",
).read()


spike_trains_list = elephant.spike_train_generation.threshold_detection(
    memb_pot_data[0].segments[0].analogsignals[0]
)

# saving to file
with open("HL23Net_0.2.spikes", "w") as f:
    gid = 0
    for spike_train in spike_trains_list:
        spikes = spike_train.as_array()
        # for file
        for s in spikes:
            print(f"{s}\t{gid}", file=f)
        gid += 1


# for plotting here
spike_data = []
gid = 0
start_cellnum = 0

for pop, cellnum in cells.items():
    all_times = []
    all_gids = []
    for spike_train in spike_trains_list[start_cellnum : start_cellnum + cellnum]:
        spikes = spike_train.as_array()
        all_times.extend(spikes)
        all_gids.extend([gid] * len(spikes))
        gid += 1

    spike_data.append({"name": pop, "times": all_times, "ids": all_gids})

    start_cellnum += cellnum

plot_spikes(spike_data=spike_data, offset=False, rates=True)
