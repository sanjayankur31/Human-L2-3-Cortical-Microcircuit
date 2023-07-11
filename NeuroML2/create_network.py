#!/usr/bin/env python3
"""
Parse population and connectivity information exported from circuit.py to
create a NeuroML representation of the network.

File: create_network.py

Copyright 2023 Ankur Sinha
"""
import logging
import pathlib

import h5py
import lems.api as lems
import neuroml
import pyneuroml
from neuroml.utils import component_factory
from pyneuroml.lems.LEMSSimulation import LEMSSimulation
from pyneuroml.neuron.nrn_export_utils import get_segment_group_name
from pyneuroml.plot.Plot import generate_plot
from pyneuroml.pynml import reload_saved_data
from pyneuroml.pynml import run_lems_with_jneuroml_neuron
from pyneuroml.pynml import write_neuroml2_file
from pyneuroml.utils import rotate_cell


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class HL23Net(object):

    """HL23 network"""

    def __init__(self, scale=0.01):
        """Init

        :param scale: network scale

        """
        object.__init__(self)

        self.network_scale = scale

        # data dumped from the simulation
        self.cell_data = h5py.File('../L23Net/Circuit_output/cell_positions_and_rotations.h5', 'r')
        self.connectivity_data = h5py.File('../L23Net/Circuit_output/synapse_connections.h5', 'r')

        # confirmed from self.cell_data.keys()
        self.cell_types = ['HL23PV', 'HL23PYR', 'HL23SST', 'HL23VIP']
        self.pop_colors = {
            'HL23PV': "0 0 1",
            'HL23PYR': "1 0 0",
            'HL23SST': "0 1 0",
            'HL23VIP': "0.5 0.5 0.5"
        }
        self.simulation_id = "HL23Sim"
        self.lems_simulation_file = "LEMS_HL23Sim.xml"

    def create_network(self):
        # set the scale of the network
        print(f"Creating network with scale {self.network_scale}")

        self.netdoc = component_factory(neuroml.NeuroMLDocument, id="HL23Network")
        self.network = self.netdoc.add(neuroml.Network, id="HL23Network", temperature="34.0 degC",
                                       notes=f"L23 network at {self.network_scale} scale", validate=False)
        self.netdoc_file_name = f"HL23Net_{self.network_scale}.net.nml"

        # synapse types
        # LEMS component definitions will be included in simulation file later
        self.synapse_components = lems.Model()
        self.synapse_components.add(lems.Include("synapses/ProbAMPANMDA.synapse.nml"))
        self.synapse_components.add(lems.Include("synapses/ProbUDF.synapse.nml"))
        self.synapse_components_file_name = "synapse_components.xml"

        self.create_cells()
        # self.create_connections()
        self.add_stimulus()

        print(self.netdoc.summary())
        self.netdoc.validate(recursive=True)

        print(f"Writing {self.synapse_components_file_name} ")
        self.synapse_components.export_to_file(self.synapse_components_file_name)

        print(f"Writing {self.netdoc_file_name} ")
        write_neuroml2_file(self.netdoc, self.netdoc_file_name, validate=False)

    def create_cells(self):
        print("Creating cells")
        # make a directory for storing rotated cells
        # we include these cells in the network document to ensure that the network
        # document doesn't get too large
        self.temp_cell_dir = "rotated_cells"
        cellfilesdir = pathlib.Path(self.temp_cell_dir)
        cellfilesdir.mkdir(exist_ok=True)

        # keep track of cell gids for our connections later, required for scaled down
        # versions when not all cells are included
        self.cell_list = []

        # create the cell populations
        for ctype in self.cell_types:
            celldataset = self.cell_data[ctype]
            # ['gid', 'x', 'y', 'z', 'x_rot', 'y_rot', 'z_rot']
            logger.debug(f"table headers are:  {celldataset.dtype.fields.keys()}")

            nml_cell = neuroml.loaders.read_neuroml2_file(f"{ctype}.cell.nml").cells[0]

            # include the cell to ensure the ion channel files are included
            self.netdoc.add(neuroml.IncludeType, href=f"{ctype}.cell.nml")

            i = 0
            step = int(1 / self.network_scale)
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
                write_neuroml2_file(rotated_cell_doc, f"{self.temp_cell_dir}/{rotated_cell.id}.cell.nml", validate=False)
                self.netdoc.add(neuroml.IncludeType, href=f"{self.temp_cell_dir}/{rotated_cell.id}.cell.nml")

                pop = self.network.add(neuroml.Population, id=f"{ctype}_pop_{gid}", type="populationList",
                                       component=rotated_cell.id)
                pop.add(neuroml.Property(tag="color", value=self.pop_colors[ctype]))
                pop.add(neuroml.Property(tag="region", value="L23"))

                pop.add(neuroml.Instance, id=0,
                        location=neuroml.Location(x=x, y=y, z=z))
                self.cell_list.append(gid)

        print(self.netdoc.summary())
        self.netdoc.validate(recursive=True)

    def create_connections(self):
        print("Creating connections")
        # count how many connections we have in total
        conn_count = 0

        # create connections
        for pretype in self.cell_types:
            for posttype in self.cell_types:
                conndataset = self.connectivity_data[f"{pretype}:{posttype}"]
                # string
                mechanism = (self.connectivity_data[f'synparams/{pretype}:{posttype}']['mechanism'][()].decode('utf-8'))

                # all ints/floats
                if "UDF" in mechanism:
                    tau_r = (self.connectivity_data[f'synparams/{pretype}:{posttype}']['tau_r'][()])
                    tau_d = (self.connectivity_data[f'synparams/{pretype}:{posttype}']['tau_d'][()])
                elif "AMPANMDA" in mechanism:
                    tau_r_AMPA = (self.connectivity_data[f'synparams/{pretype}:{posttype}']['tau_r_AMPA'][()])
                    tau_r_NMDA = (self.connectivity_data[f'synparams/{pretype}:{posttype}']['tau_r_NMDA'][()])
                    tau_d_AMPA = (self.connectivity_data[f'synparams/{pretype}:{posttype}']['tau_d_AMPA'][()])
                    tau_d_NMDA = (self.connectivity_data[f'synparams/{pretype}:{posttype}']['tau_d_NMDA'][()])
                else:
                    raise ValueError(f"Unknown mechanism found: {mechanism}")

                # common to both synapses
                Use = (self.connectivity_data[f'synparams/{pretype}:{posttype}']['Use'][()])
                Dep = (self.connectivity_data[f'synparams/{pretype}:{posttype}']['Dep'][()])
                Fac = (self.connectivity_data[f'synparams/{pretype}:{posttype}']['Fac'][()])
                gbase = (self.connectivity_data[f'synparams/{pretype}:{posttype}']['gmax'][()])
                u0 = (self.connectivity_data[f'synparams/{pretype}:{posttype}']['u0'][()])
                erev = (self.connectivity_data[f'synparams/{pretype}:{posttype}']['e'][()])

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
                                         mg="1 mM",
                                         weight_factor_NMDA="1"
                                         )

                self.synapse_components.add(syn)

                anml_cell = neuroml.loaders.read_neuroml2_file(f"{posttype}.cell.nml").cells[0]  # type: neuroml.Cell
                syn_count = 0
                cur_precell = None
                cur_postcell = None
                print(f"Creating connections: {pretype} -> {posttype} (~{int(conndataset.shape[0] * self.network_scale * self.network_scale)} conns).")

                for conn in conndataset:
                    precell = conn[0]
                    postcell = conn[1]
                    weight = conn[2]
                    delay = conn[3]
                    section = conn[4]
                    sectionx = conn[5]

                    # if both cells are not in our population, skip this connection
                    if precell not in self.cell_list or postcell not in self.cell_list:
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
                        proj = self.network.add(neuroml.Projection, id=f"proj_{precell}_{postcell}",
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

    def add_stimulus(self):
        """Add stimulus to cells.

        Currently simply adds a pulse generator input to each cell for testing
        """
        print("Adding stimuli")
        pg = self.netdoc.add(
            "PulseGenerator",
            id="pulseGen_0", delay="200ms", duration="1000ms",
            amplitude="200 pA"
        )
        ctr = 0
        for apop in self.network.populations:
            inputlist = self.network.add(
                neuroml.InputList, id=f"i{ctr}", component=pg.id, populations=apop.id,
                validate=False
            )
            inputlist.add(
                neuroml.Input, id=ctr, target=f"../{apop.id}/0/{apop.component}",
                destination="synapses"
            )
            ctr += 1

    def create_simulation(self, dt=0.025, seed=123):
        # Create simulation, and record data
        simulation = LEMSSimulation(
            sim_id=self.simulation_id, duration=2000, dt=dt,
            simulation_seed=seed
        )
        # simulation.assign_simulation_target(network.id)
        simulation.assign_simulation_target("HL23Network")
        simulation.include_neuroml2_file(f"HL23Net_{self.network_scale}.net.nml")
        # simulation.include_lems_file(self.synapse_components_file_name)
        simulation.include_lems_file("synapse_components.xml")

        simulation.create_output_file("output1", f"HL23Net_{self.network_scale}.v.dat")
        for apop in self.network.populations:
            simulation.add_column_to_output_file("output1", f"{apop.id}", f"{apop.id}/0/{apop.component}/0/v")

        simulation.save_to_file(self.lems_simulation_file)
        print(f"Saved simulation to {self.lems_simulation_file}")

    def run_sim(self):
        """Run the sim"""
        print(f"Running simulation: {self.lems_simulation_file}")
        run_lems_with_jneuroml_neuron(self.lems_simulation_file,
                                      max_memory="8G", nogui=True,
                                      show_plot_already=False)

    def plot_v_graphs(self):
        """Plot membrane potential graphs"""
        data = reload_saved_data(self.lems_simulation_file)
        logger.debug(data.keys())
        xvals = [data['t']]
        yvals = list(data.values())[1:]

        logger.debug(len(xvals * len(yvals)))
        logger.debug(yvals)

        labels = [a.split("/")[0] for a in data.keys()][1:]

        generate_plot(xvalues=xvals * len(yvals), yvalues=yvals,
                      title="Membrane potentials", labels=labels,
                      xaxis="time (ms)", yaxis="v (mV)",
                      cols_in_legend_box=2,
                      save_figure_to=f"{self.lems_simulation_file.replace('.xml', '')}_v.png")


if __name__ == "__main__":
    model = HL23Net(scale=0.01)
    model.create_network()
    model.create_simulation()
    model.run_sim()
    model.plot_v_graphs()
