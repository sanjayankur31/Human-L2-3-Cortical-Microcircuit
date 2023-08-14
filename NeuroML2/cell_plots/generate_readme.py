#!/usr/bin/env python3
"""
Generate a readme with figures
- one column per cell
- one row per property: morphology, and each ion channel


File: generate_readme.py

Copyright 2023 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""


import pathlib


cellnames = ["HL23PV", "HL23PYR", "HL23SST", "HL23VIP"]
file = "Readme.md"

pngs = [p.name for p in pathlib.Path(".").glob("*.ion.png")]

all_channels = []

for p in pngs:
    ion = p.split(".")[0].split("_", maxsplit=1)[1]
    all_channels.append(ion)

all_channels = sorted(list(set(all_channels)))
print(f"Ion channels found: {all_channels}")


with open(file, 'w') as r:
    print("# Ion channel: summary", file=r)
    print("", file=r)
    print("", file=r)

    header = "| Ion channel |"
    header_underline = "| --- |"
    for cell in cellnames:
        header += f" {cell} |"
        header_underline += " --- |"

    # print table
    print(f"{header}", file=r)
    print(f"{header_underline}", file=r)
    for c in all_channels:
        row = f"| {c} |"
        row2 = f"| {c} (dist) |"
        for cell in cellnames:
            p_name = f"{cell}_{c}.ion.png"
            if p_name in pngs:
                row += f" ![{p_name}]({p_name}) |"
                vs_dist = p_name.split(".")[0] + "_ion_vs_dist.png"
                row2 += f" ![{vs_dist}]({vs_dist}) |"
            else:
                row += " |"
                row2 += " |"
        print(row, file=r)
        if "vs_dist" in row2:
            print(row2, file=r)
