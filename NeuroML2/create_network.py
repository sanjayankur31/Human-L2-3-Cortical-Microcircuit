#!/usr/bin/env python3
"""
Parse population and connectivity information exported from circuit.py to
create a NeuroML representation of the network.

File: create_network.py

Copyright 2023 Ankur Sinha
"""

# https://docs.h5py.org/en/stable/quick.html#quick


from pyneuroml.lems.LEMSSimulation import LEMSSimulation
import pathlib
import logging

import neuroml
from neuroml.utils import component_factory
import pyneuroml
from pyneuroml.utils import rotate_cell
from pyneuroml.pynml import write_neuroml2_file
from pyneuroml.neuron.nrn_export_utils import get_segment_group_name
import lems.api as lems
import h5py


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# set the scale of the network
network_scale = 0.02
print(f"Creating network with scale {network_scale}")

cell_data = h5py.File('../L23Net/Circuit_output/cell_positions_and_rotations.h5', 'r')
# confirmed from cell_data.keys()
cell_types = ['HL23PV', 'HL23PYR', 'HL23SST', 'HL23VIP']
# cell_types = ['HL23PV']
pop_colors = {
    'HL23PV': "0 0 1",
    'HL23PYR': "1 0 0",
    'HL23SST': "0 1 0",
    'HL23VIP': "0.5 0.5 0.5"
}

netdoc = component_factory(neuroml.NeuroMLDocument, id="L23Network")
network = netdoc.add(neuroml.Network, id="L23Network", temperature="34.0 degC",
                     notes=f"L23 network at {network_scale} scale", validate=False)
netdoc_file_name = f"L23Net_{network_scale}.net.nml"

# synapse types
# LEMS component definitions will be included in simulation file later
synapse_components = lems.Model()
synapse_components.add(lems.Include("synapses/ProbAMPANMDA.synapse.nml"))
synapse_components.add(lems.Include("synapses/ProbUDF.synapse.nml"))
synapse_components_file_name = "synapse_components.xml"

# make a directory for storing rotated cells
# we include these cells in the network document to ensure that the network
# document doesn't get too large
temp_cell_dir = "rotated_cells"
cellfilesdir = pathlib.Path(temp_cell_dir)
cellfilesdir.mkdir(exist_ok=True)

# keep track of cell gids for our connections later, required for scaled down
# versions when not all cells are included
cell_list = []

# create the cell populations
for ctype in cell_types:
    celldataset = cell_data[ctype]
    # ['gid', 'x', 'y', 'z', 'x_rot', 'y_rot', 'z_rot']
    logger.debug(f"table headers are:  {celldataset.dtype.fields.keys()}")

    nml_cell = neuroml.loaders.read_neuroml2_file(f"{ctype}.cell.nml").cells[0]

    # include the cell to ensure the ion channel files are included
    netdoc.add(neuroml.IncludeType, href=f"{ctype}.cell.nml")

    i = 0
    step = int(1 / network_scale)
    for i in range(0, len(celldataset), step):
        acell = celldataset[i]
        gid = acell[0]
        x = acell[1]
        y = acell[2]
        z = acell[3]
        xrot = acell[4]
        yrot = acell[5]
        zrot = acell[6]

        rotated_cell = None
        rotated_cell = rotate_cell(nml_cell, xrot, yrot, zrot, order="xyz", relative_to_soma=True)
        rotated_cell.id = rotated_cell.id + f"_{gid}"
        rotated_cell_doc = component_factory(neuroml.NeuroMLDocument, id=f"{rotated_cell.id}_doc")
        rotated_cell_doc.add(rotated_cell)
        write_neuroml2_file(rotated_cell_doc, f"{temp_cell_dir}/{rotated_cell.id}.cell.nml", validate=False)
        netdoc.add(neuroml.IncludeType, href=f"{temp_cell_dir}/{rotated_cell.id}.cell.nml")

        pop = network.add(neuroml.Population, id=f"{ctype}_pop_{gid}", type="populationList",
                          component=rotated_cell.id)
        pop.add(neuroml.Property(tag="color", value=pop_colors[ctype]))
        pop.add(neuroml.Property(tag="region", value="L23"))

        pop.add(neuroml.Instance, id=0,
                location=neuroml.Location(x=x, y=y, z=z))
        cell_list.append(gid)

print(netdoc.summary())
netdoc.validate(recursive=True)

# count how many connections we have in total
conn_count = 0

# create connections
connectivity_data = h5py.File('../L23Net/Circuit_output/synapse_connections.h5', 'r')
for pretype in cell_types:
    for posttype in cell_types:
        conndataset = connectivity_data[f"{pretype}:{posttype}"]
        # string
        mechanism = (connectivity_data[f'synparams/{pretype}:{posttype}']['mechanism'][()].decode('utf-8'))

        # all ints/floats
        if "UDF" in mechanism:
            tau_r = (connectivity_data[f'synparams/{pretype}:{posttype}']['tau_r'][()])
            tau_d = (connectivity_data[f'synparams/{pretype}:{posttype}']['tau_d'][()])
        elif "AMPANMDA" in mechanism:
            tau_r_AMPA = (connectivity_data[f'synparams/{pretype}:{posttype}']['tau_r_AMPA'][()])
            tau_r_NMDA = (connectivity_data[f'synparams/{pretype}:{posttype}']['tau_r_NMDA'][()])
            tau_d_AMPA = (connectivity_data[f'synparams/{pretype}:{posttype}']['tau_d_AMPA'][()])
            tau_d_NMDA = (connectivity_data[f'synparams/{pretype}:{posttype}']['tau_d_NMDA'][()])
        else:
            raise ValueError(f"Unknown mechanism found: {mechanism}")

        # common to both synapses
        Use = (connectivity_data[f'synparams/{pretype}:{posttype}']['Use'][()])
        Dep = (connectivity_data[f'synparams/{pretype}:{posttype}']['Dep'][()])
        Fac = (connectivity_data[f'synparams/{pretype}:{posttype}']['Fac'][()])
        gbase = (connectivity_data[f'synparams/{pretype}:{posttype}']['gmax'][()])
        u0 = (connectivity_data[f'synparams/{pretype}:{posttype}']['u0'][()])
        erev = (connectivity_data[f'synparams/{pretype}:{posttype}']['e'][()])

        print(f"Creating synapse component: {pretype} -> {posttype}: {pretype}_{posttype}_{mechanism}.")
        if "UDF" in mechanism:
            syn = lems.Component(id_=f"{pretype}_{posttype}_{mechanism}",
                                 type_=f"{mechanism}",
                                 tau_r=f"{tau_r} ms",
                                 tau_d=f"{tau_d} ms",
                                 Use=Use,
                                 Dep=f"{Dep} ms",
                                 Fac=f"{Fac} ms",
                                 gbase=f"{gbase} uS",
                                 u0=u0,
                                 erev=f"{erev} mV",
                                 )
        elif "AMPANMDA" in mechanism:
            syn = lems.Component(id_=f"{pretype}_{posttype}_{mechanism}",
                                 type_=f"{mechanism}",
                                 tau_r_AMPA=f"{tau_r_AMPA} ms",
                                 tau_d_AMPA=f"{tau_d_AMPA} ms",
                                 tau_r_NMDA=f"{tau_r_NMDA} ms",
                                 tau_d_NMDA=f"{tau_d_NMDA} ms",
                                 Use=Use,
                                 Dep=f"{Dep} ms",
                                 Fac=f"{Fac} ms",
                                 gbase=f"{gbase} uS",
                                 u0=u0,
                                 erev=f"{erev} mV",
                                 mg="1 mM"
                                 weight_factor_NMDA="1"
                                 )

        synapse_components.add(syn)

        anml_cell = neuroml.loaders.read_neuroml2_file(f"{posttype}.cell.nml").cells[0]  # type: neuroml.Cell
        syn_count = 0
        cur_precell = None
        cur_postcell = None
        print(f"Creating connections: {pretype} -> {posttype} (~{int(conndataset.shape[0] * network_scale * network_scale)} conns).")

        for conn in conndataset:
            precell = conn[0]
            postcell = conn[1]
            weight = conn[2]
            delay = conn[3]
            section = conn[4]
            sectionx = conn[5]

            # if both cells are not in our population, skip this connection
            if precell not in cell_list or postcell not in cell_list:
                logger.debug(f"{precell} or {postcell} are not included in the network. Skipping")
                continue

            section = (section.decode("utf-8")).split(".")[1]
            neuroml_seggrp_id = get_segment_group_name(section)
            [ord_segs, cumul_lengths] = anml_cell.get_ordered_segments_in_groups(group_list=[neuroml_seggrp_id], include_cumulative_lengths=True)

            list_ord_segs = ord_segs[neuroml_seggrp_id]
            list_cumul_lengths = cumul_lengths[neuroml_seggrp_id]
            total_len = list_cumul_lengths[-1]

            section_loc = total_len * sectionx
            ind = 0
            for leng in list_cumul_lengths:
                if leng > section_loc:
                    break
                ind += 1

            post_seg = list_ord_segs[ind]

            # if it's the first segment, with ind 0, [ind - 1] is not its
            # parent segment
            if ind != 0:
                frac_along = ((section_loc - list_cumul_lengths[ind - 1]) / (list_cumul_lengths[ind] - list_cumul_lengths[ind - 1]))
                logger.debug(f"frac_along: ({section_loc} - {list_cumul_lengths[ind - 1]}) / ({list_cumul_lengths[ind]} - {list_cumul_lengths[ind - 1]}) = {frac_along}")
            else:
                frac_along = (section_loc / list_cumul_lengths[ind])
                logger.debug(f"frac_along: ({section_loc} / {list_cumul_lengths[ind]}) = {frac_along}")

            # for zero length segments
            if frac_along == -float("inf"):
                frac_along = 1

            conn_count += 1
            logger.debug(f"{conn_count}: {pretype}:{precell} -> {posttype}:{postcell} {neuroml_seggrp_id}: segment {post_seg.id} ({list_cumul_lengths[ind-1]} - {list_cumul_lengths[ind]}) at {frac_along} with mechanism {mechanism}")

            # a new projection is only required when the pre or post cell
            # change
            if precell != cur_precell or postcell != cur_postcell:
                proj = network.add(neuroml.Projection, id=f"proj_{precell}_{postcell}",
                                   presynaptic_population=f"{pretype}_pop_{precell}",
                                   postsynaptic_population=f"{posttype}_pop_{postcell}",
                                   synapse=f"{pretype}_{posttype}_{mechanism}")

                cur_precell = precell
                cur_postcell = postcell
                syn_count = 0

            try:
                proj.add(neuroml.ConnectionWD, id=syn_count,
                         pre_cell_id=f"../{pretype}_pop_{precell}/0/{pretype}_{precell}",
                         pre_segment_id=0,
                         post_cell_id=f"../{posttype}_pop_{postcell}/0/{posttype}_{postcell}",
                         post_segment_id=post_seg.id,
                         post_fraction_along=frac_along,
                         weight=weight,
                         delay=f"{delay} ms"
                         )
            except ValueError as e:
                print(f"list of cumulative lengths: {list_cumul_lengths}")
                print(f"frac_along: ({section_loc} - {list_cumul_lengths[ind - 1]}) / ({list_cumul_lengths[ind]} - {list_cumul_lengths[ind - 1]}) = {frac_along}")
                raise e

            syn_count += 1

print(netdoc.summary())
netdoc.validate(recursive=True)

print(f"Writing {synapse_components_file_name} ")
synapse_components.export_to_file(synapse_components_file_name)

print(f"Writing {netdoc_file_name} ")
write_neuroml2_file(netdoc, netdoc_file_name, validate=False)

# Create simulation, and record data
simulation_id = "HL23Sim"
simulation = LEMSSimulation(
    sim_id=simulation_id, duration=1000, dt=0.1, simulation_seed=123
)
# simulation.assign_simulation_target(network.id)
simulation.assign_simulation_target("L23Network")
simulation.include_neuroml2_file("L23Net_0.02.net.nml")
# simulation.include_lems_file(synapse_components_file_name)
simulation.include_lems_file("synapse_components.xml")

lems_simulation_file = simulation.save_to_file()
