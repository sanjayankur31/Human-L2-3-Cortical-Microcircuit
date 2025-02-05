#!/usr/bin/env python3
"""
Main script for creating scaled versions of NeuroML network and simulating it.

This parses population and connectivity information exported from circuit.py to
create a NeuroML representation of the network.

File: nml_main.py

Copyright 2023 Ankur Sinha
"""

import argparse
import bisect
import copy
import inspect
import logging
import math
import pathlib
import random
import textwrap
import time
import typing

import h5py
import lems.api as lems
import neuroml
import numpy
import pandas
import progressbar
from neuroml.loaders import read_neuroml2_file
from neuroml.utils import component_factory
from neuroml.writers import NeuroMLHdf5Writer
from pyneuroml.lems.LEMSSimulation import LEMSSimulation
from pyneuroml.neuron.nrn_export_utils import get_segment_group_name
from pyneuroml.nsgr import run_on_nsg
from pyneuroml.plot.Plot import generate_plot
from pyneuroml.plot.PlotMorphology import plot_2D
from pyneuroml.plot.PlotMorphologyVispy import plot_interactive_3D
from pyneuroml.pynml import (
    generate_sim_scripts_in_folder,
    reload_saved_data,
    run_lems_with,
    write_neuroml2_file,
)
from pyneuroml.utils import rotate_cell
from pyneuroml.utils.units import convert_to_units, get_value_in_si

logger = logging.getLogger("HL23-NeuroML2")
logger.setLevel(logging.INFO)


# https://stackoverflow.com/questions/3410976/how-to-round-a-number-to-significant-figures-in-python
# currently unused
def round_to_sig(x):
    one_less = round(x, int(math.floor(math.log10(abs(x))))) / 10
    return round(x + one_less, int(math.floor(math.log10(abs(one_less)))))


class HL23Net(object):
    """HL23 network"""

    def __init__(
        self,
        scale: float = 0.01,
        connections: bool = True,
        network_input: str = "background",
        stimulus: bool = True,
        biophysics: bool = True,
        tonic_inhibition: bool = True,
        new_cells: bool = True,
        max_memory: str = "8G",
        hdf5=True,
        rotate_cells=False,
    ):
        """Init

        :param scale: network scale
        :type scale: float
        :param connections: toggle creation of network connections
        :type connections: bool
        :param network_input: select input to provide to cells
            - "background": for OU background input
            - "step": for a constant step current
            - "none": for no input
        :type network_input: str
        :param stimulus: toggle addition of stimulus to cells for stim protocol
        :type stimulus: bool
        :param biophysics: toggle addition of biophysics to cells (otherwise
            the cells only include morphology)
        :type biophysics: bool
        :param tonic_inhibition: toggle addition of tonic inhibition to cells
            note: requires biophysics to be enabled
        :type tonic_inhibition: bool
        :param new_cells: toggle regeneration of rotated cells, otherwise
            existing rotated cells are used
        :type new_cells: bool
        :param max_memory: max memory for JVM when running pyneuroml
        :type max_memory: str
        :param hdf5: toggle exporting to HDF5 format
        :type hdf5: bool
        :param rotate_cells: generate rotated cells, useful for visualization
            but does not affect the spiking simulation (will affect LPF, but
            we're not generating those here)
        :type rotate_cells: bool

        """
        object.__init__(self)

        self.network_scale = scale
        self.connections = connections
        self.network_input = network_input
        self.stimulus = stimulus
        self.biophysics = biophysics
        self.tonic_inhibition = tonic_inhibition
        self.new_cells = new_cells
        self.max_memory = max_memory
        self.hdf5 = hdf5
        self.rotate_cells = rotate_cells
        # true by default
        self.create = True

        # data dumped from the simulation
        try:
            self.cell_data = h5py.File(
                "../L23Net/Circuit_output/cell_positions_and_rotations.h5", "r"
            )
            self.connectivity_data = h5py.File(
                "../L23Net/Circuit_output/synapse_connections.h5", "r"
            )
            self.circuit_params = pandas.read_excel(
                "../L23Net/Circuit_param.xls", sheet_name=None, index_col=0
            )
            # default synaptic parameters: are read from exported H5 files for
            # creation of synapses (above)
            self.circuit_params["syn_params"] = {
                "none": {
                    "tau_r_AMPA": 0,
                    "tau_d_AMPA": 0,
                    "tau_r_NMDA": 0,
                    "tau_d_NMDA": 0,
                    "e": 0,
                    "Dep": 0,
                    "Fac": 0,
                    "Use": 0,
                    "u0": 0,
                    "gmax": 0,
                }
            }

            # confirmed from self.cell_data.keys()
            self.cell_types = [i for i in self.circuit_params["conn_probs"].axes[0]]
        except FileNotFoundError:
            print(
                "Original simulation data files not found. Will not be able to re-create model."
            )
            self.create = False

        self.pop_colors = {
            "HL23PV": "0 0 1",
            "HL23PYR": "1 0 0",
            "HL23SST": "0 1 0",
            "HL23VIP": "0.5 0.5 0.5",
        }
        self.simulation_id = f"HL23Sim_{self.network_scale}"
        if self.rotate_cells is True:
            self.lems_simulation_file = (
                f"LEMS_HL23_{self.network_scale}_Sim.rotated.xml"
            )
        else:
            self.lems_simulation_file = f"LEMS_HL23_{self.network_scale}_Sim.xml"
        self.netdoc = None
        self.network_id = f"HL23Network_{str(self.network_scale).replace('.', '_')}"
        if self.rotate_cells is True:
            self.netdoc_file_name = f"HL23Net_{self.network_scale}.rotated.net.nml"
        else:
            self.netdoc_file_name = f"HL23Net_{self.network_scale}.net.nml"
        if self.hdf5 is True:
            self.netdoc_file_name += ".h5"
        self.lems_components_file_name = f"lems_components_{self.network_scale}.xml"
        self.stim_start = "200ms"
        self.sim_length = "1000ms"
        self.sim_end = (
            str(
                convert_to_units(
                    f"{get_value_in_si(self.stim_start) + get_value_in_si(self.sim_length)} s",
                    "ms",
                )
            )
            + "ms"
        )
        self.dt = "0.025ms"
        self.seed = 4587

    def create_network(self):
        # set the scale of the network
        if self.create is False:
            print("Not creating network")
            return

        start = time.time()
        print(f"Creating network with scale {self.network_scale}")

        self.netdoc = component_factory(neuroml.NeuroMLDocument, id="HL23Network")
        self.network = self.netdoc.add(
            neuroml.Network,
            id=self.network_id,
            temperature="34.0 degC",
            notes=f"L23 network at {self.network_scale} scale",
            validate=False,
        )

        # synapse types
        # LEMS component definitions will be included in simulation file later
        self.lems_components = lems.Model()
        self.lems_components.add(lems.Include("CaDynamics_E2_NML2.nml"))
        self.lems_components.add(lems.Include("synapses/ProbAMPANMDA.synapse.nml"))
        self.lems_components.add(lems.Include("synapses/ProbUDF.synapse.nml"))

        # Always include, set conductance to 0 if not required
        self.lems_components.add(lems.Include("channels/Tonic.nml"))

        # add all the channel definitions to
        channel_files = pathlib.Path("channels").glob("**/*.channel.nml")
        for afile in channel_files:
            logger.debug(f"Including {afile}")
            self.lems_components.add(lems.Include(str(afile)))

        self.create_cells()
        if self.connections is True:
            self.create_connections()
        if self.network_input == "background":
            self.add_background_input()
        elif self.network_input == "step":
            self.add_step_current()
        elif self.network_input == "none":
            pass
        else:
            print(f"Invalid network_input value: {self.network_input}, not adding any")
            pass
        if self.stimulus is True:
            self.add_stimulus()

        print(self.netdoc.summary())
        self.netdoc.validate(recursive=True)

        print(f"Writing {self.lems_components_file_name} ")
        self.lems_components.export_to_file(self.lems_components_file_name)

        print(f"Writing {self.netdoc_file_name} ")
        if self.hdf5:
            NeuroMLHdf5Writer.write(
                self.netdoc, self.netdoc_file_name, embed_xml=False, compress=True
            )
        else:
            write_neuroml2_file(self.netdoc, self.netdoc_file_name, validate=False)

        end = time.time()
        print(f"Creating network took: {(end - start)} seconds.")

    def create_cells(self):
        """Create cells and add them to the network.

        Each rotated cell is a new component in NeuroML, and since each
        population can only have one component attached to it when cells are
        rotated, each cell will belong to a different single celled population.

        However, for simulations where one isn't looking at LFPs etc., one does
        not need to rotate the cells. So, one cell component will be used to
        create a single population of each cell type.
        """
        start = time.time()
        print("Creating cells")
        # make a directory for storing rotated cells
        # we include these cells in the network document to ensure that the network
        # document doesn't get too large
        self.temp_cell_dir = "rotated_cells"
        cellfilesdir = pathlib.Path(self.temp_cell_dir)
        cellfilesdir.mkdir(exist_ok=True)

        # keep track of cell gids for our connections later, required for scaled down
        # versions when not all cells are included
        self.nml_cell = {}  # type: typing.Dict[str, neuroml.Cell]
        self.cell_list = {}  # type: typing.Dict[str, int]
        self.cell_list_by_type = {}
        for ctype in self.cell_types:
            self.cell_list_by_type[ctype] = []

        # create the cell populations
        for ctype in self.cell_types:
            celldataset = self.cell_data[ctype]
            # ['gid', 'x', 'y', 'z', 'x_rot', 'y_rot', 'z_rot']
            logger.debug(f"table headers are:  {celldataset.dtype.fields.keys()}")

            self.nml_cell[ctype] = neuroml.loaders.read_neuroml2_file(
                f"{ctype}.cell.nml"
            ).cells[0]  # type: neuroml.Cell
            # replace biophys with empty object
            if self.biophysics is False:
                self.nml_cell[
                    ctype
                ].biophysical_properties = neuroml.BiophysicalProperties(id="biophys")
                self.nml_cell[ctype].biophysical_properties.add(
                    neuroml.MembraneProperties, validate=False
                )
                self.nml_cell[ctype].biophysical_properties.add(
                    neuroml.IntracellularProperties
                )

            # include the cell to ensure the ion channel files are included
            # self.netdoc.add(neuroml.IncludeType, href=f"{ctype}.cell.nml")

            # If we don't rotate cells, we only need one population per cell
            # type.
            # Make a copy anyway because we're modifying the original cell
            if self.rotate_cells is False:
                unrotated_cell = copy.deepcopy(self.nml_cell[ctype])
                unrotated_cell.id = unrotated_cell.id + "_sim"
                unrotated_cell_file = (
                    f"{self.temp_cell_dir}/{unrotated_cell.id}.cell.nml"
                )
                unrotated_cell_doc = component_factory(
                    neuroml.NeuroMLDocument, id=f"{unrotated_cell.id}_doc"
                )
                unrotated_cell_doc.add(unrotated_cell)
                if self.biophysics is True:
                    self.__add_tonic_inhibition(
                        ctype, unrotated_cell, unrotated_cell_doc
                    )
                write_neuroml2_file(
                    unrotated_cell_doc, unrotated_cell_file, validate=False
                )

                self.netdoc.add(neuroml.IncludeType, href=unrotated_cell_file)

                pop = self.network.add(
                    neuroml.Population,
                    id=f"{ctype}_pop",
                    type="populationList",
                    component=unrotated_cell.id,
                )
                pop.add(neuroml.Property(tag="color", value=self.pop_colors[ctype]))
                pop.add(neuroml.Property(tag="region", value="L23"))

            i = 0
            step = int(1 / self.network_scale)
            maxcell = len(celldataset)
            # put in minimum 2 cells of each type
            if step >= maxcell:
                step = 1
                maxcell = 2
            for j in range(0, maxcell, step):
                acell = celldataset[j]
                gid = acell[0]
                x = acell[1]
                y = acell[2]
                z = acell[3]
                xrot = acell[4]
                yrot = acell[5]
                zrot = acell[6]

                if self.rotate_cells is True:
                    rotated_cell = None
                    rotated_cell_file = (
                        f"{self.temp_cell_dir}/{self.nml_cell[ctype].id}_{gid}.cell.nml"
                    )
                    rotated_cell_id = self.nml_cell[ctype].id + f"_{gid}"
                    if (
                        self.new_cells is False
                        and pathlib.Path(rotated_cell_file).exists()
                    ):
                        print(f"{rotated_cell_file} already exists, not overwriting.")
                        print("Set new_cells=True to regenerate all rotated cell files")
                    else:
                        rotated_cell = rotate_cell(
                            self.nml_cell[ctype],
                            xrot,
                            yrot,
                            zrot,
                            order="xyz",
                            relative_to_soma=True,
                        )
                        rotated_cell.id = rotated_cell_id
                        rotated_cell_doc = component_factory(
                            neuroml.NeuroMLDocument, id=f"{rotated_cell_id}_doc"
                        )

                        if self.biophysics is True:
                            self.__add_tonic_inhibition(
                                ctype, rotated_cell, rotated_cell_doc
                            )

                        rotated_cell_doc.add(rotated_cell)
                        write_neuroml2_file(
                            rotated_cell_doc, rotated_cell_file, validate=False
                        )

                    self.netdoc.add(neuroml.IncludeType, href=rotated_cell_file)

                    pop = self.network.add(
                        neuroml.Population,
                        id=f"{ctype}_pop_{gid}",
                        type="populationList",
                        component=rotated_cell_id,
                    )
                    pop.add(neuroml.Property(tag="color", value=self.pop_colors[ctype]))
                    pop.add(neuroml.Property(tag="region", value="L23"))

                    pop.add(
                        neuroml.Instance, id=0, location=neuroml.Location(x=x, y=y, z=z)
                    )
                    # track what gids we're including in our network
                    # used later when creating connections
                    self.cell_list[gid] = 0
                else:
                    pop.add(
                        neuroml.Instance, id=i, location=neuroml.Location(x=x, y=y, z=z)
                    )
                    self.cell_list[gid] = i

                i += 1

                # currently unused
                self.cell_list_by_type[ctype].append(gid)

        print(self.netdoc.summary())
        self.netdoc.validate(recursive=True)
        end = time.time()
        print(f"Creating cells took: {(end - start)} seconds")

    def __add_tonic_inhibition(self, ctype, cell, cell_doc):
        """Add tonic inhibition to provided cell

        :param ctype: cell type
        :type ctype: str
        :param cell: cell object to add to
        :type cell: neuroml.Cell
        :param cell_doc: cell doc object
        :type cell_doc: neuroml.NeuroMLDocument
        """
        if self.tonic_inhibition is True:
            cond_density = "0.000938 S_per_cm2"
        else:
            # add but deactivate
            cond_density = "0.0 S_per_cm2"

        # add gaba tonic inhibition to cells
        # definition file is included in the network, so we don't
        # re-include it here.
        if ctype == "HL23PYR" or ctype == "HL23SST":
            cell.add_channel_density(
                nml_cell_doc=cell_doc,
                cd_id="TonicInhibition",
                ion_channel="TonicPavlov2009",
                cond_density=cond_density,
                erev="-75 mV",
                group_id="all_minus_myelin",
                ion="non_specific",
                ion_chan_def_file="",
            )
        elif ctype == "HL23PV" or ctype == "HL23VIP":
            cell.add_channel_density(
                nml_cell_doc=cell_doc,
                cd_id="TonicInhibition",
                ion_channel="TonicPavlov2009",
                cond_density=cond_density,
                erev="-75 mV",
                group_id="all",
                ion="non_specific",
                ion_chan_def_file="",
            )

    def create_connections(self):
        start = time.time()
        print("Creating connections")
        # count how many connections we have in total
        conn_count = 0

        # create dicts to hold values for caching: prevents us from having to
        # get this info again and again
        ordered_segments = {}
        cumulative_lengths = {}
        for atype in self.cell_types:
            ordered_segments[f"{atype}"] = {}
            cumulative_lengths[f"{atype}"] = {}

        # create connections
        for pretype in self.cell_types:
            for posttype in self.cell_types:
                conndataset = self.connectivity_data[f"{pretype}:{posttype}"]
                # string
                mechanism = self.connectivity_data[f"synparams/{pretype}:{posttype}"][
                    "mechanism"
                ][()].decode("utf-8")

                # all ints/floats
                if "UDF" in mechanism:
                    tau_r = self.connectivity_data[f"synparams/{pretype}:{posttype}"][
                        "tau_r"
                    ][()]
                    tau_d = self.connectivity_data[f"synparams/{pretype}:{posttype}"][
                        "tau_d"
                    ][()]
                elif "AMPANMDA" in mechanism:
                    tau_r_AMPA = self.connectivity_data[
                        f"synparams/{pretype}:{posttype}"
                    ]["tau_r_AMPA"][()]
                    tau_r_NMDA = self.connectivity_data[
                        f"synparams/{pretype}:{posttype}"
                    ]["tau_r_NMDA"][()]
                    tau_d_AMPA = self.connectivity_data[
                        f"synparams/{pretype}:{posttype}"
                    ]["tau_d_AMPA"][()]
                    tau_d_NMDA = self.connectivity_data[
                        f"synparams/{pretype}:{posttype}"
                    ]["tau_d_NMDA"][()]
                else:
                    raise ValueError(f"Unknown mechanism found: {mechanism}")

                # common to both synapses
                Use = self.connectivity_data[f"synparams/{pretype}:{posttype}"]["Use"][
                    ()
                ]
                Dep = self.connectivity_data[f"synparams/{pretype}:{posttype}"]["Dep"][
                    ()
                ]
                Fac = self.connectivity_data[f"synparams/{pretype}:{posttype}"]["Fac"][
                    ()
                ]
                gbase = self.connectivity_data[f"synparams/{pretype}:{posttype}"][
                    "gmax"
                ][()]
                u0 = self.connectivity_data[f"synparams/{pretype}:{posttype}"]["u0"][()]
                erev = self.connectivity_data[f"synparams/{pretype}:{posttype}"]["e"][
                    ()
                ]

                print(
                    f"Creating synapse component: {pretype} -> {posttype}: {pretype}_{posttype}_{mechanism}."
                )
                if "UDF" in mechanism:
                    syn = lems.Component(
                        id_=f"{pretype}_{posttype}_{mechanism}",
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
                    syn = lems.Component(
                        id_=f"{pretype}_{posttype}_{mechanism}",
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
                        mg="1 mM",
                        weight_factor_NMDA="1",
                    )

                self.lems_components.add(syn)

                anml_cell = self.nml_cell[f"{posttype}"]  # type: neuroml.Cell
                syn_count = 0
                cur_precell = None
                cur_postcell = None
                print(
                    f"Creating connections: {pretype} -> {posttype} (~{int(conndataset.shape[0])} with scale of {self.network_scale})."
                )
                # if we're not rotating cells, we only need one project between
                # the single populations of different cell types
                if self.rotate_cells is False:
                    proj = component_factory(
                        neuroml.Projection,
                        id=f"proj_{pretype}_{posttype}",
                        presynaptic_population=f"{pretype}_pop",
                        postsynaptic_population=f"{posttype}_pop",
                        synapse=f"{pretype}_{posttype}_{mechanism}",
                    )  # type: neuroml.Projection

                # show a progress bar so we have some idea of what's going on
                max_value_possible = int(conndataset.shape[0])
                # rounded = round_to_sig(approx_max_value)
                bar = progressbar.ProgressBar(
                    max_val=max_value_possible, poll_interval=30
                )
                bar_count = 0
                for conn in conndataset:
                    precell = conn[0]
                    postcell = conn[1]
                    weight = conn[2]
                    delay = conn[3]
                    section = conn[4]
                    sectionx = conn[5]

                    # if both cells are not in our population, skip this connection
                    if (
                        precell not in self.cell_list.keys()
                        or postcell not in self.cell_list.keys()
                    ):
                        logger.debug(
                            f"{precell} or {postcell} are not included in the network. Skipping"
                        )
                        continue

                    bar.update(bar_count)
                    bar_count += 1

                    section = (section.decode("utf-8")).split(".")[1]
                    neuroml_seggrp_id = get_segment_group_name(section)

                    # use info if it's cached, otherwise get new info and cache
                    # it
                    try:
                        ord_segs = ordered_segments[f"{posttype}"][neuroml_seggrp_id]
                        cumul_lengths = cumulative_lengths[f"{posttype}"][
                            neuroml_seggrp_id
                        ]
                    except KeyError:
                        [
                            ord_segs,
                            cumul_lengths,
                        ] = anml_cell.get_ordered_segments_in_groups(
                            group_list=[neuroml_seggrp_id],
                            include_cumulative_lengths=True,
                        )
                        ordered_segments[f"{posttype}"][neuroml_seggrp_id] = ord_segs
                        cumulative_lengths[f"{posttype}"][neuroml_seggrp_id] = (
                            cumul_lengths
                        )

                    list_ord_segs = ord_segs[neuroml_seggrp_id]
                    list_cumul_lengths = cumul_lengths[neuroml_seggrp_id]
                    total_len = list_cumul_lengths[-1]

                    section_loc = total_len * sectionx
                    ind = bisect.bisect_left(list_cumul_lengths, section_loc)

                    post_seg = list_ord_segs[ind]

                    # if it's the first segment, with ind 0, [ind - 1] is not its
                    # parent segment
                    if ind != 0:
                        frac_along = (section_loc - list_cumul_lengths[ind - 1]) / (
                            list_cumul_lengths[ind] - list_cumul_lengths[ind - 1]
                        )
                        logger.debug(
                            f"frac_along: ({section_loc} - {list_cumul_lengths[ind - 1]}) / ({list_cumul_lengths[ind]} - {list_cumul_lengths[ind - 1]}) = {frac_along}"
                        )
                    else:
                        frac_along = section_loc / list_cumul_lengths[ind]
                        logger.debug(
                            f"frac_along: ({section_loc} / {list_cumul_lengths[ind]}) = {frac_along}"
                        )

                    # for zero length segments
                    if frac_along == -float("inf"):
                        frac_along = 1

                    conn_count += 1
                    logger.debug(
                        f"{conn_count}: {pretype}:{precell} -> {posttype}:{postcell} {neuroml_seggrp_id}: segment {post_seg.id} ({list_cumul_lengths[ind - 1]} - {list_cumul_lengths[ind]}) at {frac_along} with mechanism {mechanism}"
                    )

                    # a new projection is only required when the pre or post cell
                    # change
                    if self.rotate_cells is True:
                        if precell != cur_precell or postcell != cur_postcell:
                            proj = self.network.add(
                                neuroml.Projection,
                                id=f"proj_{precell}_{postcell}",
                                presynaptic_population=f"{pretype}_pop_{precell}",
                                postsynaptic_population=f"{posttype}_pop_{postcell}",
                                synapse=f"{pretype}_{posttype}_{mechanism}",
                            )

                            cur_precell = precell
                            cur_postcell = postcell
                            syn_count = 0

                        try:
                            proj.add(
                                neuroml.ConnectionWD,
                                id=syn_count,
                                pre_cell_id=f"../{pretype}_pop_{precell}/0/{pretype}_{precell}",
                                pre_segment_id=0,
                                post_cell_id=f"../{posttype}_pop_{postcell}/0/{posttype}_{postcell}",
                                post_segment_id=post_seg.id,
                                post_fraction_along=frac_along,
                                weight=str(float(weight) / self.network_scale),
                                delay=f"{delay} ms",
                            )
                        except ValueError as e:
                            print(f"list of cumulative lengths: {list_cumul_lengths}")
                            print(
                                f"frac_along: ({section_loc} - {list_cumul_lengths[ind - 1]}) / ({list_cumul_lengths[ind]} - {list_cumul_lengths[ind - 1]}) = {frac_along}"
                            )
                            raise e
                    else:
                        try:
                            proj.add(
                                neuroml.ConnectionWD,
                                id=syn_count,
                                pre_cell_id=f"../{pretype}_pop/{self.cell_list[precell]}/{pretype}_sim",
                                pre_segment_id=0,
                                post_cell_id=f"../{posttype}_pop/{self.cell_list[postcell]}/{posttype}_sim",
                                post_segment_id=post_seg.id,
                                post_fraction_along=frac_along,
                                weight=str(float(weight) / self.network_scale),
                                delay=f"{delay} ms",
                            )
                        except ValueError as e:
                            print(f"list of cumulative lengths: {list_cumul_lengths}")
                            print(
                                f"frac_along: ({section_loc} - {list_cumul_lengths[ind - 1]}) / ({list_cumul_lengths[ind]} - {list_cumul_lengths[ind - 1]}) = {frac_along}"
                            )
                            raise e

                    syn_count += 1
                bar.finish()

                # Only add projection to network if there's at least one
                # connection in it.
                # Not required with rotated cells because projection creation
                # and connection creation go hand in hand
                if self.rotate_cells is False:
                    if len(proj.connection_wds) > 0:
                        self.network.add(proj)

        end = time.time()
        print(f"Creating connections took: {(end - start)} seconds.")

    def add_background_input(self):
        """Add background input to cells."""
        start = time.time()
        # Component Type definition
        self.lems_components.add(lems.Include("synapses/Gfluct.nml"))

        # base excitatory conductances (from paper, given here in uS)
        cell_type_ge0 = {
            "HL23PYR": 0.000028,
            "HL23PV": 0.00028,
            "HL23SST": 0.00003,
            "HL23VIP": 0.000066,
        }
        # store input locations per cell type and create input components for
        # each cell type also, since all cells are indentical
        # values: { 'rel distance' : [(seg id, frac along)] }
        cell_type_input_locations = {}  # type: typing.Dict[str, typing.Dict[float, typing.Tuple[int, float]]]
        for cell_type, cell in self.nml_cell.items():
            cell_type_input_locations[cell_type] = {}
            extremeties = cell.get_extremeties()
            try:
                basal_segs = cell.get_all_segments_in_group("basal_dendrite_group")
            except Exception:
                # PV cell doesn't have basal, only a dendrite group
                basal_segs = cell.get_all_segments_in_group("dendrite_group")

            # input points on basal dendrites
            longest_basal_branch_length = 0
            for segid, distance in extremeties.items():
                if distance > longest_basal_branch_length and segid in basal_segs:
                    longest_basal_branch_length = distance

            half_way_basal = longest_basal_branch_length / 2
            segs_basal = cell.get_segments_at_distance(half_way_basal)

            # create input component for 0.5
            g_e0 = cell_type_ge0[cell_type] * math.exp(0.5)
            gfluct_component = lems.Component(
                id_=f"Gfluct_{cell_type}_basal_0_5",
                type_="Gfluct",
                start=self.stim_start,
                stop=self.sim_end,
                dt=self.dt,
                E_e="0mV",
                E_i="-80mV",
                g_e0=f"{g_e0} uS",
                g_i0="0pS",
                tau_e="65ms",
                tau_i="20ms",
                std_e=f"{g_e0} uS",
                std_i="0pS",
            )
            self.lems_components.add(gfluct_component)

            # get segments to place input at
            cell_type_input_locations[cell_type][0.5] = []
            for seg, frac_along in segs_basal.items():
                if seg in basal_segs:
                    cell_type_input_locations[cell_type][0.5].append((seg, frac_along))

        # create a new input list for each population, and each location
        # because input list takes a component as an argument, and a different
        # component is required for each location
        input_list_ctr = 0
        input_ctr = 0
        for pop in self.network.populations:
            # cell name
            cell_type = pop.component.split("_")[0]
            input_segs = cell_type_input_locations[cell_type]

            for rel_dist, seginfos in input_segs.items():
                # one input list per population per component
                inputlist = self.network.add(
                    "InputList",
                    id=f"Gfluct_{input_list_ctr}",
                    component=f"Gfluct_{cell_type}_basal_{str(rel_dist).replace('.', '_')}",
                    populations=pop.id,
                    validate=False,
                )
                input_list_ctr += 1
                for seg, frac_along in seginfos:
                    if self.rotate_cells is True:
                        inputlist.add(
                            "Input",
                            id=f"{input_ctr}",
                            target=f"../{pop.id}/0/{pop.component}",
                            destination="synapses",
                            segment_id=seg,
                            fraction_along=frac_along,
                        )
                        input_ctr += 1
                    else:
                        for inst in pop.instances:
                            inputlist.add(
                                "Input",
                                id=f"{input_ctr}",
                                target=f"../{pop.id}/{inst.id}/{pop.component}",
                                destination="synapses",
                                segment_id=seg,
                                fraction_along=frac_along,
                            )
                            input_ctr += 1

        # additional for pyr apical
        cell_type = "HL23PYR"
        cell_type_input_locations = {}  # type: typing.Dict[str, typing.Dict[float, typing.Tuple[int, float]]]
        cell_type_input_locations[cell_type] = {}
        pyr_cell = self.nml_cell[cell_type]
        extremeties = pyr_cell.get_extremeties()
        longest_apical_branch_length = 0
        pyr_apical_segs = pyr_cell.get_all_segments_in_group("apical_dendrite_group")
        for segid, distance in extremeties.items():
            if distance > longest_apical_branch_length and segid in pyr_apical_segs:
                longest_apical_branch_length = distance

        apical_input_distances = [0.1, 0.3, 0.5, 0.7, 0.9]
        for d in apical_input_distances:
            # create the input component:
            g_e0 = cell_type_ge0[cell_type] * math.exp(d)
            # create input component for use at each distance point
            gfluct_component = lems.Component(
                id_=f"Gfluct_HL23PYR_apical_{str(d).replace('.', '_')}",
                type_="Gfluct",
                start=self.stim_start,
                stop=self.sim_end,
                dt=self.dt,
                E_e="0mV",
                E_i="-80mV",
                g_e0=f"{g_e0} uS",
                g_i0="0pS",
                tau_e="65ms",
                tau_i="20ms",
                std_e=f"{g_e0} uS",
                std_i="0pS",
            )
            self.lems_components.add(gfluct_component)

            # get segments to place at
            cell_type_input_locations[cell_type][d] = []
            segs_apical = pyr_cell.get_segments_at_distance(
                d * longest_apical_branch_length
            )
            for seg, frac_along in segs_apical.items():
                if seg in pyr_apical_segs:
                    cell_type_input_locations[cell_type][d].append((seg, frac_along))

        # create inputs and input lists for pyr apical
        input_segs = cell_type_input_locations[cell_type]

        for pop in self.network.populations:
            if not "PYR" in pop.id:
                continue
            for rel_dist, seginfos in input_segs.items():
                # one input list per population per component
                inputlist = self.network.add(
                    "InputList",
                    id=f"Gfluct_apical_{input_list_ctr}",
                    component=f"Gfluct_{cell_type}_apical_{str(rel_dist).replace('.', '_')}",
                    populations=pop.id,
                    validate=False,
                )
                input_list_ctr += 1
                for seg, frac_along in seginfos:
                    if self.rotate_cells is True:
                        inputlist.add(
                            "Input",
                            id=f"{input_ctr}",
                            target=f"../{pop.id}/0/{pop.component}",
                            destination="synapses",
                            segment_id=seg,
                            fraction_along=frac_along,
                        )
                        input_ctr += 1
                    else:
                        for inst in pop.instances:
                            inputlist.add(
                                "Input",
                                id=f"{input_ctr}",
                                target=f"../{pop.id}/{inst.id}/{pop.component}",
                                destination="synapses",
                                segment_id=seg,
                                fraction_along=frac_along,
                            )
                            input_ctr += 1

        end = time.time()
        print(f"Adding background input took: {(end - start)} seconds.")

    # TODO: incomplete, to be done after bg input implementation
    def add_stimulus(self):
        """Add additional stimulus to cells."""

        logger.setLevel(logging.DEBUG)
        print("Adding stimuli")
        # cell_nums = [self.circuit_params['SING_CELL_PARAM'].at['cell_num', name] for name in self.cell_types]
        # from circuit.py
        input_ctr = 0
        for acell_type in self.cell_types:
            # iterate over rows
            for row in self.circuit_params["STIM_PARAM"].axes[0]:
                # initialise parameters
                num_cells = 0  # number of cells to provide input to
                start_index = 0  # unused, is set to 0
                num_stim = 0  # number of spikes
                interval = 0  # spike interval
                start_time = 0  # spikes start time
                delay = 0  # spikes delay
                delay_range = 0  # spikes delay range
                loc_num = 0  # number of locations on cell to provide spike to
                loc = "all"  # what locations to provide spikes to
                gmax = 0.0  # gmax
                stim_type = ""  # type of synapse for stimulus
                syn_params = ""  # parameters of synapse

                for col in self.circuit_params["STIM_PARAM"].axes[1]:
                    # found the cell row
                    if (
                        "cell_name" == col
                        and self.circuit_params["STIM_PARAM"].at[row, col] == acell_type
                    ):
                        num_cells = int(
                            self.circuit_params["STIM_PARAM"].at[row, "num_cells"]
                            * self.network_scale
                        )
                        start_index = self.circuit_params["STIM_PARAM"].at[
                            row, "start_index"
                        ]
                        num_stim = self.circuit_params["STIM_PARAM"].at[row, "num_stim"]
                        interval = self.circuit_params["STIM_PARAM"].at[row, "interval"]
                        # currently unused: NeuroML does not have a
                        # SpikeGenerator that allows setting a start time
                        start_time = self.circuit_params["STIM_PARAM"].at[
                            row, "start_time"
                        ]
                        delay = self.circuit_params["STIM_PARAM"].at[row, "delay"]
                        delay_range = self.circuit_params["STIM_PARAM"].at[
                            row, "delay_range"
                        ]
                        loc_num = self.circuit_params["STIM_PARAM"].at[row, "loc_num"]
                        # loc: unused: "dend", which corresponds to all
                        # dendritic segments
                        loc = self.circuit_params["STIM_PARAM"].at[row, "loc"]
                        gmax = self.circuit_params["STIM_PARAM"].at[row, "gmax"]
                        stim_type = self.circuit_params["STIM_PARAM"].at[
                            row, "stim_type"
                        ]
                        syn_params = self.circuit_params["STIM_PARAM"].at[
                            row, "syn_params"
                        ]

                        # load the single template cell, since choosing sections does
                        # not depend on rotation
                        a_nml_cell = self.nml_cell[acell_type]  # type: neuroml.Cell

                        # loc is always "dend"
                        logger.debug(f"loc: {loc}")
                        dendritic_segments_ids = a_nml_cell.get_all_segments_in_group(
                            "dendrite_group"
                        )

                        # dendritic_segments = [nml_cell.get_segment(seg) for seg in dendritic_segments_ids]
                        segment_areas = numpy.array(
                            [
                                a_nml_cell.get_segment_surface_area(seg)
                                for seg in dendritic_segments_ids
                            ]
                        )
                        segment_areas = segment_areas / segment_areas.sum()

                        cells = self.cell_list_by_type[acell_type]
                        ctr = 0
                        # if it's too small a model, at least give one cell input
                        if num_cells == 0:
                            num_cells = 1
                        while ctr <= num_cells:
                            cell_id = f"{acell_type}_{cells[ctr]}"
                            pop_id = f"{acell_type}_pop_{cells[ctr]}"

                            # get segments
                            input_segments = numpy.random.choice(
                                dendritic_segments_ids,
                                size=loc_num,
                                replace=False,
                                p=segment_areas,
                            )
                            for aseg in input_segments:
                                print(
                                    f"Adding sg_{input_ctr}: i_{input_ctr}/{input_ctr} to {pop_id}/{cell_id}/{aseg}"
                                )
                                # currently SpikeGenerator does not allow a start time, so
                                # this is unused
                                time_delay = 0
                                while time_delay <= 0:
                                    time_delay = numpy.random.uniform(
                                        low=delay, high=delay + delay_range
                                    )
                                gen = self.netdoc.add(
                                    neuroml.SpikeGenerator,
                                    id=f"sg_{input_ctr}",
                                    period=f"{interval} ms",
                                )
                                input_list = self.network.add(
                                    neuroml.InputList,
                                    id=f"i_{input_ctr}",
                                    component=gen.id,
                                    populations=pop_id,
                                )
                                input_list.add(
                                    neuroml.Input,
                                    id=f"{input_ctr}",
                                    target=f"../{pop_id}/0/{cell_id}",
                                    destination="synapses",
                                    segment_id=aseg,
                                )
                                input_ctr += 1

                            ctr += 1

                        break
        logger.setLevel(logging.INFO)

    def visualize_network(self, min_cells=5):
        """Generate morph plots"""
        print(
            f"Visualising network at scale {self.network_scale} with {min_cells} in full"
        )
        # if the network has been constructed, copy over the populations and
        # include the cells
        if self.netdoc is not None:
            nml_file = neuroml.NeuroMLDocument(id="network_copy")
            nml_file.add(neuroml.Network, id="dummy_network")
            # copy over populations
            nml_file.networks[0].populations = self.netdoc.networks[0].populations
            nml_file.includes = self.netdoc.includes
            # include cells
            for inc in self.nml_file.includes:
                incfile = read_neuroml2_file(inc.href)
                for cell in incfile.cells:
                    nml_file.add(cell)
        else:
            # otherwise read the file in once so it's not read in repeatedly
            print(f"Reading {self.netdoc_file_name}")
            nml_file = read_neuroml2_file(self.netdoc_file_name)

        PYR_cells = []
        SST_cells = []
        PV_cells = []
        VIP_cells = []
        PYR_point_cells = []
        SST_point_cells = []
        PV_point_cells = []
        VIP_point_cells = []

        # because each cell is in a separate population
        print("Picking point cells")
        for inc in nml_file.includes:
            cell = inc.href.split("/")[1].replace(".cell.nml", "")
            if "PYR" in cell:
                PYR_cells.append(cell)
            if "PV" in cell:
                PV_cells.append(cell)
            if "VIP" in cell:
                VIP_cells.append(cell)
            if "SST" in cell:
                SST_cells.append(cell)

        print(f"Got {len(PYR_cells)} PYR cells")
        print(f"Got {len(SST_cells)} SST cells")
        print(f"Got {len(PV_cells)} PV cells")
        print(f"Got {len(VIP_cells)} VIP cells")

        # at least plot min_cells cells of each type, and all the rest as points
        if (len(PYR_cells) - min_cells) > 0:
            PYR_point_cells = random.sample(PYR_cells, len(PYR_cells) - min_cells)
        if (len(SST_cells) - min_cells) > 0:
            SST_point_cells = random.sample(SST_cells, len(SST_cells) - min_cells)
        if (len(PV_cells) - min_cells) > 0:
            PV_point_cells = random.sample(PV_cells, len(PV_cells) - min_cells)
        if (len(VIP_cells) - min_cells) > 0:
            VIP_point_cells = random.sample(VIP_cells, len(VIP_cells) - min_cells)
        print(
            f"Picked {len(PYR_point_cells) + len(SST_point_cells) + len(PV_point_cells) + len(VIP_point_cells)} point cells"
        )

        """
        print("Removing axons from PYR cells")
        for ac in nml_file.cells:
            if ac.id not in PYR_point_cells:
                axons = ac.get_all_segments_in_group(ac.get_segment_group("axon_group"))
                for s in ac.morphology.segments:
                    if s.id in axons:
                        ac.morphology.segments.remove(s)

        """
        """
        for plane in ["xy", "yz", "zx"]:
            print(f"Plotting {plane} with {point_fraction} fraction as point cells")
            plot_2D(
                nml_file=nml_file,
                plane2d=plane,
                min_width=1,
                nogui=True,
                title="",
                plot_type="constant",
                save_to_file=f"{self.netdoc_file_name.replace('.nml', '')}.{plane}.png",
                plot_spec={
                    "point_cells": PYR_point_cells + SST_point_cells + VIP_point_cells + PV_point_cells
                }
            )

            print(f"Plotting {plane} with all cells as point cells")
            plot_2D(
                nml_file=nml_file,
                plane2d=plane,
                min_width=1,
                nogui=True,
                title="",
                plot_type="constant",
                save_to_file=f"{self.netdoc_file_name.replace('.nml', '')}.{plane}.allpoints.png",
                plot_spec={
                    "point_fraction": 1.0
                }
            )

            print(f"Plotting {plane} with a single cells of each type in detail")
            plot_2D(
                nml_file=nml_file,
                plane2d=plane,
                min_width=1,
                nogui=True,
                title="",
                plot_type="constant",
                save_to_file=f"{self.netdoc_file_name.replace('.nml', '')}.{plane}.points.png",
                plot_spec={
                    "point_cells": PYR_cells[1:] + SST_cells[1:] + VIP_cells[1:] + PV_cells[1:]
                }
            )
            """
        plot_interactive_3D(
            self.netdoc_file_name,
            plot_type="detailed",
            plot_spec={
                "point_cells": PYR_point_cells
                + SST_point_cells
                + VIP_point_cells
                + PV_point_cells
            },
            min_width=1.0,
        )

    def create_simulation(self, dt=None, seed=123, record_data=True):
        """Create simulation, record data.

        :param dt: override dt value (in ms)
        :type dt: str
        :param seed: seed
        :type seed: int
        :param record_data: toggle whether data should be recorded
        :type record_data: bool
        """
        start = time.time()
        if dt is None:
            dt = self.dt

        if seed is None:
            seed = self.seed

        simulation = LEMSSimulation(
            sim_id=self.simulation_id,
            duration=float(self.sim_length.replace("ms", "")),
            dt=float(dt.replace("ms", "")),
            simulation_seed=seed,
        )
        simulation.assign_simulation_target(self.network_id)
        network_file = pathlib.Path(f"HL23Net_{self.network_scale}.net.nml")
        if not network_file.is_file():
            logger.error(
                f"Network file not found at HL23Net_{self.network_scale}.net.nml."
            )
            logger.error("Please use the --create_network flag to create a network")
            return

        self.network_doc = read_neuroml2_file(
            f"HL23Net_{self.network_scale}.net.nml", include_includes=False
        )
        self.network = self.network_doc.networks[0]

        simulation.include_neuroml2_file(f"HL23Net_{self.network_scale}.net.nml")
        simulation.include_lems_file(self.lems_components_file_name)

        if record_data is True:
            simulation.create_output_file(
                "output1", f"HL23Net_{self.network_scale}.v.dat"
            )
            print(f"Saving data to: HL23Net_{self.network_scale}.v.dat")
            for apop in self.network.populations:
                for inst in apop.instances:
                    simulation.add_column_to_output_file(
                        "output1",
                        f"{apop.id}_{inst.id}",
                        f"{apop.id}/{inst.id}/{apop.component}/0/v",
                    )

        simulation.save_to_file(self.lems_simulation_file)
        print(f"Saved simulation to {self.lems_simulation_file}")
        end = time.time()
        print(f"Creating simulation took: {(end - start)} seconds.")

    def run_sim(
        self,
        engine: str = "jneuroml_neuron",
        cluster: typing.Union[str, bool] = False,
        **kwargs,
    ):
        """Run the sim

        :param engine: engine to use (jneuroml_neuron/jneuroml_netpyne)
        :type engine: str
        :param cluster: toggle submitting to nsg or qsub
            Use "nsg_dry" for a dry NSG run (won't submit)
            Use "qsub" to generate a qsub script instead
        :type cluster: bool or str
        :param **kwargs: other engine + nsg specific args
        """
        if cluster is False:
            print(f"Running simulation: {self.lems_simulation_file}")
            run_lems_with(
                engine,
                self.lems_simulation_file,
                max_memory=self.max_memory,
                nogui=True,
                show_plot_already=False,
                **kwargs,
            )
        elif cluster == "qsub":
            print(f"Generating files but not running: {self.lems_simulation_file}")
            # TODO: create in a different folder like NSG
            tdir = generate_sim_scripts_in_folder(
                engine=engine,
                lems_file_name=self.lems_simulation_file,
                root_dir=".",
                max_memory=self.max_memory,
                nogui=True,
                show_plot_already=False,
                **kwargs,
            )
            # remove trailing backslash
            if tdir[-1] == "/":
                tdir = tdir[:-1]
            qsub_fn = f"{tdir}/{tdir}_generated/{self.lems_simulation_file}.sh"
            netpyne_simfile = self.lems_simulation_file.replace(".xml", "_netpyne.py")
            print(f"Generating qsub script for use on cluster: {qsub_fn}")
            with open(qsub_fn, "w") as f:
                print(
                    inspect.cleandoc(
                        textwrap.dedent(
                            f"""
                            #!/bin/bash -l
                            #$ -pe mpi 36
                            #$ -l mem=4G
                            #$ -l h_rt=6:00:00
                            #$ -cwd
                            #$ -m be

                            source ~/.venv311/bin/activate
                            nrnivmodl
                            gerun python3 {netpyne_simfile} -nogui
                            """
                        )
                    ),
                    file=f,
                )
        elif cluster == "nsg_dry":
            print(
                f"Preparing to run on NSG (but not submitting): {self.lems_simulation_file}"
            )
            run_on_nsg(
                engine,
                self.lems_simulation_file,
                dry_run=True,
                max_memory=self.max_memory,
                **kwargs,
            )
        elif cluster == "nsg":
            print(f"Running simulation on NSG: {self.lems_simulation_file}")
            run_on_nsg(
                engine, self.lems_simulation_file, max_memory=self.max_memory, **kwargs
            )
            print("Press Ctrl C to interrupt the waiting process")
            print(
                "You can use `nsg_job -l` to check the status of the job and download the simulation reults"
            )

    def plot_v_graphs(self):
        """Plot membrane potential graphs"""
        data = reload_saved_data(
            self.lems_simulation_file,
            plot=True,
            show_plot_already=True,
            reload_events=False,
        )
        """
        logger.debug(data.keys())
        xvals = [data["t"]]
        yvals = list(data.values())[1:]

        logger.debug(len(xvals * len(yvals)))
        logger.debug(yvals)

        labels = [a.split("/")[0] for a in data.keys()][1:]

        generate_plot(
            xvalues=xvals * len(yvals),
            yvalues=yvals,
            title="Membrane potentials",
            labels=labels,
            xaxis="time (ms)",
            yaxis="v (mV)",
            cols_in_legend_box=2,
            save_figure_to=f"{self.lems_simulation_file.replace('.xml', '')}_v.png",
        )
        """

    def add_step_current(self):
        """Add a constant step current to all cells.

        Useful for testing the cells to ensure they produce correct behaviour,
        in an unconnected network, for example

        """
        print("Adding step current to each cell")
        # a single pulse generator
        pg = self.netdoc.add(
            neuroml.PulseGenerator,
            id="pg_0",
            delay="100ms",
            duration=self.sim_length,
            amplitude="0.2nA",
        )
        input_ctr = 0
        for pop in self.network.populations:
            input_list = self.network.add(
                "InputList", id=f"input_{pop.id}", component=pg.id, populations=pop.id
            )
            for inst in pop.instances:
                input_list.add(
                    "Input",
                    id=input_ctr,
                    target=f"../{pop.id}/{inst.id}/{pop.component}",
                    destination="synapses",
                )
                input_ctr += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="NeuroML_YaoEtAl2022",
        description="Generate and simulate NeuroML version of Yao et al 2022",
    )
    # general options
    general_args = parser.add_argument_group("General options")
    general_args.add_argument(
        "-s",
        "--scale",
        action="store",
        help="Scale of network",
        default="1",
        required=True,
        type=float,
    )
    general_args.add_argument(
        "-n",
        "--new_cells",
        action="store_true",
        help="Create and use new cells instead of re-using existing cells",
        default=False,
    )
    general_args.add_argument(
        "-b",
        "--biophysics",
        action="store_false",
        help="Disable inclusion of biophysics",
        default=True,
    )
    general_args.add_argument(
        "-t",
        "--tonic_inhibition",
        action="store_true",
        help="Enable tonic inhibition in the network",
        default=False,
    )
    general_args.add_argument(
        "-c",
        "--connections",
        action="store_false",
        help="Disable creation of network connections",
        default=True,
    )
    general_args.add_argument(
        "-i",
        "--network_input",
        action="store",
        help="Type of input to give to network",
        choices=["background", "step", "none"],
        default="background",
    )
    general_args.add_argument(
        "-l",
        "--stimulus",
        action="store_true",
        help="Enable external stimulus",
        default=False,
    )
    general_args.add_argument(
        "-5",
        "--hdf5",
        action="store_true",
        help="Enable exporting as HDF5",
        default=False,
    )
    general_args.add_argument(
        "-r",
        "--rotate_cells",
        action="store_true",
        help="Enable rotation of cells: not required for simulation since placements of cells in space does not affect it",
        default=False,
    )

    # actions
    action_args = parser.add_argument_group("Actions")
    action_args.add_argument(
        "--create_network",
        action="store_true",
        help="Create a new network",
        default=False,
    )
    action_args.add_argument(
        "--create_simulation",
        action="store_true",
        help="Create new simulation",
        default=False,
    )
    action_args.add_argument(
        "--visualize_network",
        action="store",
        metavar="<min number of cells to visualise with full morphology>",
        help="Visualise network with provided minimum number of cells",
        type=int,
        default=-1,
    )

    run_parser = action_args.add_mutually_exclusive_group()
    run_parser.add_argument(
        "--simulate_neuron",
        action="store_true",
        help="Simulation using single threaded NEURON simulation",
    )
    run_parser.add_argument(
        "--simulate_netpyne_nsg",
        action="store_true",
        help="Generate configuration to run using NetPyNE on NSG",
    )
    run_parser.add_argument(
        "--simulate_netpyne_qsub",
        action="store_true",
        help="Generate configuration to run using NetPyNE on a cluster using qsub",
    )

    # NSG options
    nsg_args = parser.add_argument_group("NSG options")

    # parse
    nsg_args.add_argument(
        "--number_cores", action="store", default="64", help="Number of cores requested"
    )
    nsg_args.add_argument(
        "--number_nodes", action="store", default="8", help="Number of nodes"
    )
    nsg_args.add_argument(
        "--tasks_per_node",
        action="store",
        default="64",
        help="Number of tasks per node (cores * number nodes?)",
    )
    nsg_args.add_argument(
        "--memory_per_node",
        action="store",
        default="240",
        help="Memory requested per node in GB",
    )
    nsg_args.add_argument(
        "--runtime", action="store", default="5", help="Run time requested in hours"
    )
    nsg_args.add_argument(
        "--environment",
        action="store",
        default="OSBv2_EXPANSE_0_7_3",
        help="Environment to request",
    )

    args = parser.parse_args()

    if args.scale < 0.01:
        logger.warning(
            "A scale < 0.01 is not recommended because the small number of neurons will result in no connections between populations at all."
        )
        logger.warning(
            "Please `pynml -graph 4d <network file>` or `pynml -matrix 1 <network file>` to see connectivity information"
        )

    model = HL23Net(
        scale=args.scale,
        new_cells=args.new_cells,
        biophysics=args.biophysics,
        tonic_inhibition=args.tonic_inhibition,
        connections=args.connections,
        network_input=args.network_input,
        stimulus=args.stimulus,
        hdf5=args.hdf5,
        rotate_cells=args.rotate_cells,
    )
    if args.create_network:
        model.create_network()
    if args.create_simulation:
        model.create_simulation()
    if args.visualize_network != -1:
        model.visualize_network(min_cells=args.visualize_network)

    # simulation
    if args.simulate_neuron:
        model.run_sim(engine="jneuroml_neuron", cluster=False, skip_run=False)
    elif args.simulate_netpyne_nsg:
        model.run_sim(
            engine="jneuroml_netpyne",
            cluster="nsg",
            nsg_sim_config={
                "number_cores_": args.number_cores,
                "tasks_per_node_": args.tasks_per_node,
                "number_nodes_": args.number_nodes,
                "number_gbmemorypernode_": args.memory_per_node,
                "runtime_": args.runtime,
                "toolId": args.environment,
                "nrnivmodl_o_": "1",
                "cmdlineopts_": "-nogui",
            },
            nogui=True,
        )
    elif args.simulate_netpyne_qsub:
        model.run_sim(engine="jneuroml_netpyne", cluster="qsub")
