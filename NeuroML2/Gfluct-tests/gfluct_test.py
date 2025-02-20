#!/usr/bin/env python3
"""
Test Gfluct inputs on single cells

File: gfluct_test.py

Copyright 2025 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""

import logging
import math
import time
from types import MappingProxyType

import lems.api as lems
import neuroml
from frozendict import frozendict
from neuroml.utils import component_factory
from pyneuroml.io import read_neuroml2_file, write_neuroml2_file
from pyneuroml.lems.LEMSSimulation import LEMSSimulation
from pyneuroml.plot.PlotMorphology import plot_2D
from pyneuroml.plot.PlotMorphologyVispy import plot_interactive_3D
from pyneuroml.plot.PlotTimeSeries import (
    plot_time_series,
    plot_time_series_from_lems_file,
)
from pyneuroml.runners import run_lems_with_jneuroml_neuron
from pyneuroml.utils.units import convert_to_units, get_value_in_si, split_nml2_quantity

logger = logging.getLogger("gfluct_test")
logger.setLevel(logging.DEBUG)


class GfluctTest(object):
    """Class for Gfluct test"""

    def __init__(self):
        """No-op"""
        self.simulation_id = "GfluctTest"
        self.lems_components = lems.Model()
        self.connections = False
        self.nml_cell = {}  # type: typing.Dict[str, neuroml.Cell]
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
        self.rotate_cells = False

    def loadcell(self):
        """Load a cell to test
        :returns: TODO

        """
        celldoc = read_neuroml2_file("./HL23PYR.cell.nml")
        cell = celldoc.cells[0]
        self.nml_cell["HL23PYR"] = cell

        self.netdoc = component_factory(neuroml.NeuroMLDocument, id="HL23Network")
        self.netdoc.add(neuroml.IncludeType, href="./HL23PYR.cell.nml")
        self.netdoc_file_name = "GfluctTest.net.nml"
        self.lems_simulation_file = "LEMS_GfluctTest.xml"
        self.lems_components_file_name = "lems_components.xml"
        self.network_id = f"{cell.id}_network"
        self.network = self.netdoc.add(
            neuroml.Network,
            id=self.network_id,
            temperature="34.0 degC",
            validate=False,
        )

        pop = self.network.add(
            neuroml.Population,
            id=f"{cell.id}_pop",
            component=cell.id,
            type="populationList",
        )
        pop.add(
            neuroml.Instance,
            id="0",
            location=component_factory(neuroml.Location, x=0, y=0, z=0),
        )

    def add_background_input(self):
        """Add background input to cells.

        Taken from main nml_main.py script. Will be updated here and added back
        there when finally ready

        Note that this component was validated here:
        https://github.com/OpenSourceBrain/StochasticityShowcase/pull/9
        """
        start = time.time()
        # Component Type definition
        self.lems_components.add(lems.Include("Gfluct.nml"))

        # base excitatory conductances (from paper, given here in uS)
        cell_type_ge0 = {
            "HL23PYR": "0.012 uS",
            "HL23PV": "280 pS",
            "HL23SST": "30 pS",
            "HL23VIP": "66 pS",
        }
        cell_type_gi0 = {
            "HL23PYR": "0.058 uS",
            "HL23PV": "140 pS",
            "HL23SST": "15 pS",
            "HL23VIP": "33 pS",
        }
        std_e = "0.00001 uS"
        std_i = "0.00001 uS"
        tau_e = "2.96 ms"
        tau_i = "9.6 ms"

        # store input locations per cell type and create input components for
        # each cell type also, since all cells are indentical
        # values: { 'rel distance' : [(seg id, frac along)] }
        self.cell_type_input_locations = {}  # type: typing.Dict[str, typing.Dict[float, typing.Tuple[int, float]]]
        for cell_type, cell in self.nml_cell.items():
            self.cell_type_input_locations[cell_type] = {}
            extremeties = cell.get_extremeties()
            logger.debug(
                f"Cell extremeties for {cell_type} cell are ({len(extremeties)}): {extremeties}"
            )
            try:
                basal_segs = cell.get_all_segments_in_group("basal_dendrite_group")
            except Exception:
                # PV cell doesn't have basal, only a dendrite group
                basal_segs = cell.get_all_segments_in_group("dendrite_group")

            # Note that the paper says: "placed at halfway the length of each
            # dendritic arbor to ensure similar levels of inputs along each
            # dendritic path to the soma". However, the code in
            # "net_functions.hoc" does not do that. It finds the length of the
            # longest branch and places the processes at sites that are at half
            # this length on all branches, not at the mid point of each branch.
            # So, we follow the code, and not the paper here.

            # input points on basal dendrites
            longest_basal_branch_length = 0
            for segid, distance in extremeties.items():
                if distance > longest_basal_branch_length and segid in basal_segs:
                    longest_basal_branch_length = distance
            logger.info(
                f"Longest basal branch for {cell_type} cell is {longest_basal_branch_length}"
            )

            half_way_basal = longest_basal_branch_length / 2
            segs_halfway = cell.get_segments_at_distance(half_way_basal)
            self.input_segs = {}

            count_inputs = 100
            picked = 0
            for sg, frac in segs_halfway.items():
                if sg in basal_segs:
                    if sg in [671, 1437, 73, 1315]:
                        continue
                    print(f"Picked basal seg: {cell.get_segment_location_info(sg)}")
                    print(
                        f"Picked basal seg: {cell.get_segment_volume(sg)}, {cell.get_segment_surface_area(sg)}"
                    )
                    self.input_segs[sg] = frac
                    picked += 1
                    if picked >= count_inputs:
                        break

            # self.input_segs = {0: 0}

            logger.debug(
                f"Basal input sites for {cell_type} cell are ({len(self.input_segs)}): {self.input_segs}"
            )

            cell_spec = frozendict(
                {
                    # 'marker_size': (5, 5),
                    "marker_size": 5,
                    "marker_color": "red",
                }
            )

            # plot the marked dendrites
            hi_spec = {cell.id: {"cell_color": "blue"}}
            for sg in self.input_segs.keys():
                hi_spec[cell.id][sg] = cell_spec

            plot_2D(cell, highlight_spec=hi_spec)
            # plot_interactive_3D(cell, highlight_spec=hi_spec)

            # create input component for 0.5
            g_e0, unite = split_nml2_quantity(cell_type_ge0[cell_type])
            g_i0, uniti = split_nml2_quantity(cell_type_gi0[cell_type])
            std_e0, unite1 = split_nml2_quantity(std_e)
            std_i0, uniti1 = split_nml2_quantity(std_i)
            gfluct_component = lems.Component(
                id_=f"Gfluct_{cell_type}_basal_0_5",
                type_="Gfluct",
                start=self.stim_start,
                stop=self.sim_end,
                dt=self.dt,
                E_e="0mV",
                E_i="-80mV",
                g_e0=f"{g_e0 * math.exp(0.5)} {unite}",
                g_i0=f"{g_i0 * math.exp(0.5)} {uniti}",
                tau_e=tau_e,
                tau_i=tau_i,
                # std_e=f"{std_e0 * math.exp(0.5)} {unite1}",
                std_e=std_e,
                # std_i=f"{std_i0 * math.exp(0.5)} {uniti1}",
                std_i=std_i,
            )
            self.lems_components.add(gfluct_component)

            # get segments to place input at
            self.cell_type_input_locations[cell_type][0.5] = []
            for seg, frac_along in self.input_segs.items():
                self.cell_type_input_locations[cell_type][0.5].append((seg, frac_along))

        # create a new input list for each population, and each location
        # because input list takes a component as an argument, and a different
        # component is required for each location
        input_list_ctr = 0
        input_ctr = 0
        for pop in self.network.populations:
            # cell name
            cell_type = pop.component.split("_")[0]
            cell_input_segs = self.cell_type_input_locations[cell_type]

            for rel_dist, seginfos in cell_input_segs.items():
                # one input list per population per component
                inputlist = self.network.add(
                    "InputList",
                    id=f"Gfluct_basal_{input_list_ctr}",
                    component=f"Gfluct_{cell_type}_basal_{str(rel_dist).replace('.', '_')}",
                    populations=pop.id,
                    validate=False,
                )
                input_list_ctr += 1
                for aseg in seginfos:
                    print(f"Doing seg {aseg}")
                    seg, frac_along = aseg
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
        """
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
            g_e0, unit = split_nml2_quantity(cell_type_ge0[cell_type])
            # create input component for use at each distance point
            gfluct_component = lems.Component(
                id_=f"Gfluct_HL23PYR_apical_{str(d).replace('.', '_')}",
                type_="Gfluct",
                start=self.stim_start,
                stop=self.sim_end,
                dt=self.dt,
                E_e="0mV",
                E_i="-80mV",
                g_e0=f"{g_e0 * math.exp(d)} {unit}",
                g_i0=g_i0,
                tau_e=tau_e,
                tau_i=tau_i,
                std_e=std_e if std_e else f"{g_e0 * math.exp(d)} {unit}",
                std_i=std_i,
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

        """
        end = time.time()
        print(f"Adding background input took: {(end - start)} seconds.")

    def write_to_file(self):
        """Write model and simulation to file"""
        print(f"Writing {self.lems_components_file_name} ")
        self.lems_components.export_to_file(self.lems_components_file_name)

        print(f"Writing {self.netdoc_file_name} ")
        write_neuroml2_file(self.netdoc, self.netdoc_file_name, validate=False)

        simulation = LEMSSimulation(
            sim_id=self.simulation_id,
            duration=float(self.sim_length.replace("ms", "")),
            dt=float(self.dt.replace("ms", "")),
            simulation_seed=self.seed,
        )
        simulation.assign_simulation_target(self.network_id)

        simulation.include_neuroml2_file(f"{self.netdoc_file_name}")
        simulation.include_lems_file(self.lems_components_file_name)

        simulation.create_output_file("output1", "Gfluct.v.dat")
        for apop in self.network.populations:
            for inst in apop.instances:
                simulation.add_column_to_output_file(
                    "output1",
                    f"{apop.id}_{inst.id}",
                    f"{apop.id}/{inst.id}/{apop.component}/0/v",
                )
                for reldist, seginfo in self.cell_type_input_locations[
                    "HL23PYR"
                ].items():
                    for aseg in seginfo:
                        seg, frac = aseg
                        simulation.add_column_to_output_file(
                            "output1",
                            f"{apop.id}_{inst.id}_{seg}",
                            f"{apop.id}/{inst.id}/{apop.component}/{seg}/v",
                        )

        simulation.save_to_file(self.lems_simulation_file)
        print(f"Saved simulation to {self.lems_simulation_file}")

        data = run_lems_with_jneuroml_neuron(
            self.lems_simulation_file,
            nogui=True,
            load_saved_data=True,
            show_plot_already=False,
        )

        plot_time_series(data, offset=False, labels=True)
        # plot_time_series_from_lems_file(self.lems_simulation_file)


if __name__ == "__main__":
    obj = GfluctTest()
    obj.loadcell()
    obj.add_background_input()
    obj.write_to_file()
