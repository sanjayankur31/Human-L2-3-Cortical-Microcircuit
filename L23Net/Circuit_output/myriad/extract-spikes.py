#!/usr/bin/env python3
"""
Extract spikes from pickled file to text files

File: L23Net/Circuit_output/myriad/extract-spikes.py

Copyright 2024 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""

import numpy

var = numpy.load("./SPIKES_Seed1234.pickle", allow_pickle=True)
dict_var = var.item()
print(dict_var.keys())

print(len(dict_var["times"][3]))
print(len(dict_var["gids"][3]))

# for each population
with open("pop-spikes.txt", "w") as f:
    spikes = {}
    for i in range(0, 4):
        times = dict_var["times"][i]
        gids = dict_var["gids"][i]

        # {spike time: gid}
        for j in range(0, len(gids)):
            gid = gids[j]
            for s in times[j]:
                spikes[s] = gid

        sorted_spike = dict(sorted(spikes.items()))

    for spike, cell in sorted_spike.items():
        print(f"{cell}\t{spike}", file=f)
