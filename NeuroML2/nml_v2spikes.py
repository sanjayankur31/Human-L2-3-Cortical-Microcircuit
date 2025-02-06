#!/usr/bin/env python3
"""
Get and plot spikes from membrane potential data, using neo and elephant.

File: NeuroML2/nml_v2spikes.py

Copyright 2025 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""

import elephant
import quantities
from neo.io import AsciiSignalIO
from pyneuroml.plot.PlotSpikes import plot_spikes

# TODO
# Should be able to get number of cells from the NeuroML network description
# instead of hard coding it.
# The order of recordings will have to be picked from the LEMS file, though
cells = {"PYR": 160, "SST": 10, "PV": 14, "VIP": 16}
# 160 PYR, 10 SST, 14 PV, 16 VIP
data_filename = "HL23Net_0.2.v.dat"
# 200 colums, but skip the 0th time column, and the spurious last column
cols = list(range(1, 201))

# sampling rate is 1/dt
memb_pot_data = AsciiSignalIO(
    "HL23Net_0.2.v.dat",
    timecolumn=None,
    sampling_rate=40000.0 * quantities.Hz,
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
