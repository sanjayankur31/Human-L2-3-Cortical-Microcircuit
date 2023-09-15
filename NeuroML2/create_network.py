#!/usr/bin/env python3
"""
Parse population and connectivity information exported from circuit.py to
create a NeuroML representation of the network.

File: create_network.py

Copyright 2023 Ankur Sinha
"""
import math
import copy
import logging
import pathlib
import sys
import typing

import h5py
import lems.api as lems
import neuroml
import numpy
import pandas
import pyneuroml
from neuroml.loaders import read_neuroml2_file
from neuroml.utils import component_factory
from pyneuroml.lems.LEMSSimulation import LEMSSimulation
from pyneuroml.neuron.nrn_export_utils import get_segment_group_name
from pyneuroml.plot.Plot import generate_plot
from pyneuroml.plot.PlotMorphology import plot_2D
from pyneuroml.pynml import (
    reload_saved_data,
    write_neuroml2_file,
    run_lems_with
)
from pyneuroml.utils import rotate_cell
from pyneuroml.nsgr import run_on_nsg

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


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

        """
        object.__init__(self)

        self.network_scale = scale
        self.connections = connections
        self.network_input = network_input
        self.stimulus = stimulus
        self.biophysics = biophysics
        self.tonic_inhibition = tonic_inhibition
        self.new_cells = new_cells

        # data dumped from the simulation
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
        self.pop_colors = {
            "HL23PV": "0 0 1",
            "HL23PYR": "1 0 0",
            "HL23SST": "0 1 0",
            "HL23VIP": "0.5 0.5 0.5",
        }
        self.simulation_id = f"HL23Sim_{self.network_scale}"
        self.lems_simulation_file = f"LEMS_HL23_{self.network_scale}_Sim.xml"
        self.netdoc = None
        self.network_id = "HL23Network"
        self.netdoc_file_name = f"HL23Net_{self.network_scale}.net.nml"
        self.lems_components_file_name = f"lems_components_{self.network_scale}.xml"
        self.sim_length = "1000ms"
        self.dt = "0.025ms"
        self.seed = 4587

    def create_network(self):
        # set the scale of the network
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
        if self.tonic_inhibition:
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
        write_neuroml2_file(self.netdoc, self.netdoc_file_name, validate=False)

    def create_cells(self):
        """Create all rotated cells and add them to the network.

        Each rotated cell is added to a separate population because in NeuroML,
        a population can only include a single Component, but each rotated cell
        is a different component.
        """
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
        self.cell_list = []
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
            ).cells[
                0
            ]  # type: neuroml.Cell
            # replace biophys with empty object
            if self.biophysics is False:
                self.nml_cell[
                    ctype
                ].biophysical_properties = neuroml.BiophysicalProperties(id="biophys")
                self.nml_cell[ctype].biophysical_properties.add(
                    neuroml.MembraneProperties
                )
                self.nml_cell[ctype].biophysical_properties.add(
                    neuroml.IntracellularProperties
                )

            # include the cell to ensure the ion channel files are included
            # self.netdoc.add(neuroml.IncludeType, href=f"{ctype}.cell.nml")

            i = 0
            step = int(1 / self.network_scale)
            maxcell = len(celldataset)
            # put in minimum 2 cells of each type
            if step >= maxcell:
                step = 1
                maxcell = 2
            for i in range(0, maxcell, step):
                acell = celldataset[i]
                gid = acell[0]
                x = acell[1]
                y = acell[2]
                z = acell[3]
                xrot = acell[4]
                yrot = acell[5]
                zrot = acell[6]

                rotated_cell = None
                rotated_cell_file = (
                    f"{self.temp_cell_dir}/{self.nml_cell[ctype].id}_{gid}.cell.nml"
                )
                rotated_cell_id = self.nml_cell[ctype].id + f"_{gid}"
                if self.new_cells is False and pathlib.Path(rotated_cell_file).exists():
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

                    if self.biophysics is True and self.tonic_inhibition is True:
                        # add gaba tonic inhibition to cells
                        # definition file is included in the network, so we don't
                        # re-include it here.
                        if ctype == "HL23PYR" or ctype == "HL23SST":
                            rotated_cell.add_channel_density(
                                nml_cell_doc=rotated_cell_doc,
                                cd_id="TonicInhibition",
                                ion_channel="TonicPavlov2009",
                                cond_density="0.000938 S_per_cm2",
                                erev="-75 mV",
                                group_id="all_minus_myelin",
                                ion="non_specific",
                                ion_chan_def_file="",
                            )
                        elif ctype == "HL23PV" or ctype == "HL23VIP":
                            rotated_cell.add_channel_density(
                                nml_cell_doc=rotated_cell_doc,
                                cd_id="TonicInhibition",
                                ion_channel="TonicPavlov2009",
                                cond_density="0.000938 S_per_cm2",
                                erev="-75 mV",
                                group_id="all",
                                ion="non_specific",
                                ion_chan_def_file="",
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
                self.cell_list.append(gid)
                self.cell_list_by_type[ctype].append(gid)

        print(self.netdoc.summary())
        self.netdoc.validate(recursive=True)

    def create_connections(self):
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
                    f"Creating connections: {pretype} -> {posttype} (~{int(conndataset.shape[0] * self.network_scale * self.network_scale)} conns)."
                )

                for conn in conndataset:
                    precell = conn[0]
                    postcell = conn[1]
                    weight = conn[2]
                    delay = conn[3]
                    section = conn[4]
                    sectionx = conn[5]

                    # if both cells are not in our population, skip this connection
                    if precell not in self.cell_list or postcell not in self.cell_list:
                        logger.debug(
                            f"{precell} or {postcell} are not included in the network. Skipping"
                        )
                        continue

                    section = (section.decode("utf-8")).split(".")[1]
                    neuroml_seggrp_id = get_segment_group_name(section)

                    # use info if it's cached, otherwise get new info and cache
                    # it
                    try:
                        ord_segs = ordered_segments[f"{posttype}"][neuroml_seggrp_id]
                        cumul_lengths = cumulative_lengths[f"{posttype}"][neuroml_seggrp_id]
                    except KeyError:
                        [
                            ord_segs,
                            cumul_lengths,
                        ] = anml_cell.get_ordered_segments_in_groups(
                            group_list=[neuroml_seggrp_id], include_cumulative_lengths=True
                        )
                        ordered_segments[f"{posttype}"][neuroml_seggrp_id] = ord_segs
                        cumulative_lengths[f"{posttype}"][neuroml_seggrp_id] = cumul_lengths

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
                        f"{conn_count}: {pretype}:{precell} -> {posttype}:{postcell} {neuroml_seggrp_id}: segment {post_seg.id} ({list_cumul_lengths[ind-1]} - {list_cumul_lengths[ind]}) at {frac_along} with mechanism {mechanism}"
                    )

                    # a new projection is only required when the pre or post cell
                    # change
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
                            weight=weight,
                            delay=f"{delay} ms",
                        )
                    except ValueError as e:
                        print(f"list of cumulative lengths: {list_cumul_lengths}")
                        print(
                            f"frac_along: ({section_loc} - {list_cumul_lengths[ind - 1]}) / ({list_cumul_lengths[ind]} - {list_cumul_lengths[ind - 1]}) = {frac_along}"
                        )
                        raise e

                    syn_count += 1

    def add_background_input(self):
        """Add background input to cells."""
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
            gfluct_component = lems.Component(id_=f"Gfluct_{cell_type}_0_5",
                                              type_="Gfluct", start="0ms",
                                              stop=self.sim_length, dt=self.dt,
                                              E_e="0mV",
                                              E_i="-80mV", g_e0=f"{g_e0} uS", g_i0="0pS",
                                              tau_e="65ms", tau_i="20ms",
                                              std_e=f"{g_e0} uS", std_i="0pS")
            self.lems_components.add(gfluct_component)

            # get segments to place input at
            cell_type_input_locations[cell_type][0.5] = []
            for seg, frac_along in segs_basal.items():
                if seg in basal_segs:
                    cell_type_input_locations[cell_type][0.5].append((seg, frac_along))

        # additional for pyr apical cells
        pyr_cell = self.nml_cell["HL23PYR"]
        extremeties = pyr_cell.get_extremeties()
        longest_apical_branch_length = 0
        pyr_apical_segs = pyr_cell.get_all_segments_in_group("apical_dendrite_group")
        for segid, distance in extremeties.items():
            if distance > longest_apical_branch_length and segid in pyr_apical_segs:
                longest_apical_branch_length = distance

        apical_input_distances = [0.1, 0.3, 0.5, 0.7, 0.9]
        for d in apical_input_distances:
            # create the input component
            g_e0 = cell_type_ge0[cell_type] * math.exp(d)
            # create input component for use at each distance point
            gfluct_component = lems.Component(id_=f"Gfluct_HL23PYR_{str(d).replace('.', '_')}",
                                              type_="Gfluct", start="0ms",
                                              stop=self.sim_length, dt=self.dt,
                                              E_e="0mV",
                                              E_i="-80mV", g_e0=f"{g_e0} uS", g_i0="0pS",
                                              tau_e="65ms", tau_i="20ms",
                                              std_e=f"{g_e0} uS", std_i="0pS")
            self.lems_components.add(gfluct_component)

            # get segments to place at
            segs_apical = pyr_cell.get_segments_at_distance(d * longest_apical_branch_length)
            for seg, frac_along in segs_apical.items():
                if seg in pyr_apical_segs:
                    try:
                        cell_type_input_locations["HL23PYR"][d].append((seg, frac_along))
                    except KeyError:
                        cell_type_input_locations["HL23PYR"][d] = []
                        cell_type_input_locations["HL23PYR"][d].append((seg, frac_along))

        # create a new input list for each population, and each location
        # because input list takes a component as an argument, and a different
        # component is required for each location
        input_list_ctr = 0
        input_ctr = 0
        for pop in self.network.populations:
            # cell name
            cell_type = pop.component.split("_")[0]
            # temporarily skip cells I haven't sorted out yet
            try:
                input_segs = cell_type_input_locations[cell_type]
            except KeyError:
                continue

            for rel_dist, seginfos in input_segs.items():
                # one input list per population per component
                inputlist = self.network.add("InputList", id=f"Gfluct_{input_list_ctr}",
                                             component=f"Gfluct_{cell_type}_{str(rel_dist).replace('.', '_')}",
                                             populations=pop.id,
                                             validate=False)
                input_list_ctr += 1
                for (seg, frac_along) in seginfos:
                    inputlist.add("Input", id=f"{input_ctr}",
                                  target=f"../{pop.id}/0/{pop.component}",
                                  destination="synapses", segment_id=seg,
                                  fraction_along=frac_along)
                    input_ctr += 1

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

    def visualize_network(self):
        """Generate morph plots"""
        # if the network has been constructed, use the network doc object, but
        # and add the included cells so we also get morphologies plotted
        if self.netdoc is not None:
            nml_file = copy.deepcopy(self.netdoc)
            for inc in nml_file.includes:
                incfile = read_neuroml2_file(inc.href)
                for cells in incfile.cells:
                    nml_file.add(cells)

            nml_file.includes = []
        else:
            # otherwise read the file in once so it's not read in repeatedly
            nml_file = read_neuroml2_file(self.netdoc_file_name, include_includes=True)

        for plane in ["xy", "yz", "zx"]:
            print(f"Plotting {plane}")
            plot_2D(
                nml_file=nml_file,
                plane2d=plane,
                min_width=4,
                nogui=True,
                title="",
                plot_type="constant",
                save_to_file=f"{self.netdoc_file_name.replace('.nml', '')}.{plane}.png",
            )

    def create_simulation(self, dt=None, seed=123):
        """Create simulation, record data"""
        if dt is None:
            dt = self.dt

        if seed is None:
            seed = self.seed

        simulation = LEMSSimulation(
            sim_id=self.simulation_id,
            duration=float(self.sim_length.replace("ms", "")),
            dt=float(dt.replace("ms", "")),
            simulation_seed=seed
        )
        simulation.assign_simulation_target(self.network_id)
        simulation.include_neuroml2_file(f"HL23Net_{self.network_scale}.net.nml")
        simulation.include_lems_file(self.lems_components_file_name)

        simulation.create_output_file("output1", f"HL23Net_{self.network_scale}.v.dat")
        for apop in self.network.populations:
            simulation.add_column_to_output_file(
                "output1", f"{apop.id}", f"{apop.id}/0/{apop.component}/0/v"
            )

        simulation.save_to_file(self.lems_simulation_file)
        print(f"Saved simulation to {self.lems_simulation_file}")

    def run_sim(self, engine: str = "jneuroml_neuron", nsg: typing.Union[str, bool] = False):
        """Run the sim"""
        if nsg is False:
            print(f"Running simulation: {self.lems_simulation_file}")
            run_lems_with(
                engine,
                self.lems_simulation_file,
                max_memory="8G",
                nogui=True,
                show_plot_already=False,
            )
        elif nsg == "dry":
            print(f"Preparing to run on NSG (but not submitting): {self.lems_simulation_file}")
            run_on_nsg(engine, self.lems_simulation_file,
                       dry_run=True, max_memory="8G")
        else:
            print(f"Running simulation on NSG: {self.lems_simulation_file}")
            run_on_nsg(engine, self.lems_simulation_file,
                       max_memory="8G")

    def plot_v_graphs(self):
        """Plot membrane potential graphs"""
        data = reload_saved_data(self.lems_simulation_file)
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

    def add_step_current(self):
        """Add a constant step current to all cells.

        Useful for testing the cells to ensure they produce correct behaviour,
        in an unconnected network, for example

        """
        print("Adding step current to each cell")
        # a single pulse generator
        pg = self.netdoc.add(neuroml.PulseGenerator, id="pg_0",
                             delay="100ms", duration=self.sim_length,
                             amplitude="0.2nA")
        for pop in self.network.populations:
            self.network.add(neuroml.ExplicitInput, target=f"{pop.id}[0]",
                             input=pg.id)


if __name__ == "__main__":
    scale = 0.02
    if len(sys.argv) == 2:
        scale = float(sys.argv[1])
    elif len(sys.argv) > 2:
        print("Only one argument accepted, using first as scale value")
        scale = float(sys.argv[1])

    model = HL23Net(
        scale=scale,
        new_cells=False,
        biophysics=True,
        tonic_inhibition=True,
        connections=True,
        network_input="background",
        stimulus=False,
    )
    model.create_network()
    # model.visualize_network()
    model.create_simulation()
    model.run_sim()
    # model.plot_v_graphs()
