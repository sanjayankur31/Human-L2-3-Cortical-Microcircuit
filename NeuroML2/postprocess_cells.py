#!/usr/bin/env python3
"""
Post process and add biophysics to cells.

We make any updates to the morphology, and add biophysics.

File: NeuroML2/postprocess_cells.py

Copyright 2022 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""


import neuroml
from neuroml.loaders import read_neuroml2_file
from neuroml.writers import NeuroMLWriter
from pyneuroml.analysis import generate_current_vs_frequency_curve


def load_and_clean_cell(cellname: str):
    """Load a cell, and clean it to prepare it for further modifications.

    :param cellname: name of cell.
        the file containing the cell should then be <cell>.morph.cell.nml
    :returns: document with cell
    :rtype: neuroml.NeuroMLDocument

    """
    celldoc = read_neuroml2_file(f"{cellname}.morph.cell.nml")  # type: neuroml.NeuroMLDocument
    celldoc.networks = []
    cell = celldoc.cells[0]  # type: neuroml.Cell
    cell.setup_nml_cell(use_convention=False)
    cell.optimise_segment_groups()
    cell.id = cellname
    cell.notes = cell.notes.replace("NeuronTemplate_0_0", cellname)
    cell.notes += ". Reference: Yao, H. K.; Guet-McCreight, A.; Mazza, F.; Moradi Chameh, H.; Prevot, T. D.; Griffiths, J. D.; Tripathy, S. J.; Valiante, T. A.; Sibille, E. & Hay, E.  Reduced inhibition in depression impairs stimulus processing in human cortical microcircuits Cell Reports, Elsevier, 2022, 38"

    # clean annotations from segment groups
    for sg in cell.morphology.segment_groups:
        sg.annotation = None

    return celldoc


# "HL23PV" "HL23PYR" "HL23SST" "HL23VIP"
def postprocess_HL23PV():
    """Post process HL23PV and add biophysics"""
    cellname = "HL23PV"
    celldoc = load_and_clean_cell(cellname)
    cell = celldoc.cells[0]  # type: neuroml.Cell

    """
    print("Segment groups:")
    for sg in cell.morphology.segment_groups:
        print(f"** {sg.id} **")
        print(cell.get_all_segments_in_group(sg.id))
        print()
    """
    # biophysics
    # include calcium dynamics component
    celldoc.add(neuroml.IncludeType(href="CaDynamics_E2_NML2.nml"), validate=False)

    # all
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="pas",
                             ion_channel="pas",
                             cond_density="0.00011830111773572024 S_per_cm2",
                             erev="-83.92924122901199 mV",
                             group_id="all",
                             ion="pas",
                             ion_chan_def_file="pas.channel.nml")
    cell.set_resistivity("0.1 kohm_cm", group_id="all")
    cell.set_specific_capacitance("2 uF_per_cm2")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Ih",
                             ion_channel="Ih",
                             cond_density="2.7671764064314368e-05 S_per_cm2",
                             erev="-45 mV",
                             group_id="all",
                             ion="hcn",
                             ion_chan_def_file="Ih.channel.nml")

    # somatic
    sgs = cell.get_segment_groups_by_substring("soma")
    for sgid, sg in sgs.items():
        # Na
        cell.add_channel_density(nml_cell_doc=celldoc,
                                 cd_id="NaTg_somatic",
                                 ion_channel="NaTg",
                                 cond_density="0.49958525078702043 S_per_cm2",
                                 erev="50 mV",
                                 group_id=sgid,
                                 ion="Na",
                                 ion_chan_def_file="NaTg.channel.nml")
        cell.add_channel_density(nml_cell_doc=celldoc,
                                 cd_id="Nap_somatic",
                                 ion_channel="Nap",
                                 cond_density="0.008795461417521086 S_per_cm2",
                                 erev="50 mV",
                                 group_id=sgid,
                                 ion="Na",
                                 ion_chan_def_file="Nap.channel.nml")

        # K
        cell.add_channel_density(nml_cell_doc=celldoc,
                                 cd_id="K_P_somatic",
                                 ion_channel="K_P",
                                 cond_density="9.606092478937705e-06 S_per_cm2",
                                 erev="-85 mV",
                                 group_id=sgid,
                                 ion="K",
                                 ion_chan_def_file="K_P.channel.nml")
        cell.add_channel_density(nml_cell_doc=celldoc,
                                 cd_id="K_T_somatic",
                                 ion_channel="K_T",
                                 cond_density="0.0011701702607527396 S_per_cm2",
                                 erev="-85 mV",
                                 group_id=sgid,
                                 ion="K",
                                 ion_chan_def_file="K_T.channel.nml")
        cell.add_channel_density(nml_cell_doc=celldoc,
                                 cd_id="Kv3_1_somatic",
                                 ion_channel="Kv3_1",
                                 cond_density="2.9921080101237565 S_per_cm2",
                                 erev="-85 mV",
                                 group_id=sgid,
                                 ion="K",
                                 ion_chan_def_file="Kv3_1.channel.nml")
        cell.add_channel_density(nml_cell_doc=celldoc,
                                 cd_id="Im_somatic",
                                 ion_channel="Im",
                                 cond_density="0.04215865946497755 S_per_cm2",
                                 erev="-85 mV",
                                 group_id=sgid,
                                 ion="K",
                                 ion_chan_def_file="Im.channel.nml")
        cell.add_channel_density(nml_cell_doc=celldoc,
                                 cd_id="SK_somatic",
                                 ion_channel="SK",
                                 cond_density="3.7265770903193036e-06 S_per_cm2",
                                 erev="-85 mV",
                                 group_id=sgid,
                                 ion="K",
                                 ion_chan_def_file="SK.channel.nml")
        # Ca
        # https://www.neuron.yale.edu/neuron/static/new_doc/modelspec/programmatic/ions.html
        cell.add_channel_density_v(
            "ChannelDensityNernst",
            nml_cell_doc=celldoc,
            id="Ca_HVA_somatic",
            ion_channel="Ca_HVA",
            cond_density="0.00017953651378188165 S_per_cm2",
            segment_groups=sgid,
            ion="ca",
            ion_chan_def_file="Ca_HVA.channel.nml")
        cell.add_channel_density_v(
            "ChannelDensityNernst",
            nml_cell_doc=celldoc,
            id="Ca_LVA_somatic",
            ion_channel="Ca_LVA",
            cond_density="0.09250008555398015 S_per_cm2",
            segment_groups=sgid,
            ion="ca",
            ion_chan_def_file="Ca_LVA.channel.nml")
        # internal and external concentrations are set to defaults that NEURON
        # starts with
        cell.add_intracellular_property("Species", validate=False,
                                        id="ca",
                                        concentration_model="CaDynamics_E2_NML2_PV_somatic",
                                        ion="ca",
                                        initial_concentration="5.0E-11 mol_per_cm3",
                                        initial_ext_concentration="2.0E-6 mol_per_cm3",
                                        segment_groups=sgid)

    # axonal
    sgs = cell.get_segment_groups_by_substring("axon")
    for sgid, sg in sgs.items():
        # Na
        cell.add_channel_density(nml_cell_doc=celldoc,
                                 cd_id="NaTg_axonal",
                                 ion_channel="NaTg",
                                 cond_density="0.10914576408883477 S_per_cm2",
                                 erev="50 mV",
                                 group_id=sgid,
                                 ion="Na",
                                 ion_chan_def_file="NaTg.channel.nml")
        cell.add_channel_density(nml_cell_doc=celldoc,
                                 cd_id="Nap_axonal",
                                 ion_channel="Nap",
                                 cond_density="0.001200899579358837 S_per_cm2",
                                 erev="50 mV",
                                 group_id=sgid,
                                 ion="Na",
                                 ion_chan_def_file="Nap.channel.nml")

        # K
        cell.add_channel_density(nml_cell_doc=celldoc,
                                 cd_id="K_P_axonal",
                                 ion_channel="K_P",
                                 cond_density="0.6854776593761795 S_per_cm2",
                                 erev="-85 mV",
                                 group_id=sgid,
                                 ion="K",
                                 ion_chan_def_file="K_P.channel.nml")
        cell.add_channel_density(nml_cell_doc=celldoc,
                                 cd_id="K_T_axonal",
                                 ion_channel="K_T",
                                 cond_density="0.07603372775662909 S_per_cm2",
                                 erev="-85 mV",
                                 group_id=sgid,
                                 ion="K",
                                 ion_chan_def_file="K_T.channel.nml")
        cell.add_channel_density(nml_cell_doc=celldoc,
                                 cd_id="Kv3_1_axonal",
                                 ion_channel="Kv3_1",
                                 cond_density="2.988867483754507 S_per_cm2",
                                 erev="-85 mV",
                                 group_id=sgid,
                                 ion="K",
                                 ion_chan_def_file="Kv3_1.channel.nml")
        cell.add_channel_density(nml_cell_doc=celldoc,
                                 cd_id="Im_axonal",
                                 ion_channel="Im",
                                 cond_density="0.029587905136596156 S_per_cm2",
                                 erev="-85 mV",
                                 group_id=sgid,
                                 ion="K",
                                 ion_chan_def_file="Im.channel.nml")
        cell.add_channel_density(nml_cell_doc=celldoc,
                                 cd_id="SK_axonal",
                                 ion_channel="SK",
                                 cond_density="0.5121938998281017 S_per_cm2",
                                 erev="-85 mV",
                                 group_id=sgid,
                                 ion="K",
                                 ion_chan_def_file="SK.channel.nml")
        # Ca
        cell.add_channel_density_v(
            "ChannelDensityNernst",
            nml_cell_doc=celldoc,
            id="Ca_HVA_somatic",
            ion_channel="Ca_HVA",
            cond_density="0.002961469262723619 S_per_cm2",
            segment_groups=sgid,
            ion="ca",
            ion_chan_def_file="Ca_HVA.channel.nml")
        cell.add_channel_density_v(
            "ChannelDensityNernst",
            nml_cell_doc=celldoc,
            id="Ca_LVA_axonal",
            ion_channel="Ca_LVA",
            cond_density="5.9457835817342756e-05 S_per_cm2",
            segment_groups=sgid,
            ion="ca",
            ion_chan_def_file="Ca_LVA.channel.nml")
        # internal and external concentrations are set to defaults that NEURON
        # starts with
        cell.add_intracellular_property("Species", validate=False,
                                        id="ca",
                                        concentration_model="CaDynamics_E2_NML2_PV_axonal",
                                        ion="ca",
                                        initial_concentration="5.0E-11 mol_per_cm3",
                                        initial_ext_concentration="2.0E-6 mol_per_cm3",
                                        segment_groups=sgid)
    cell.validate(recursive=True)
    NeuroMLWriter.write(celldoc, f"{cellname}.cell.nml", close=True)

    generate_current_vs_frequency_curve(nml2_file=f"{cellname}.cell.nml",
                                        cell_id=cellname,
                                        simulator="jNeuroML_NEURON")


if __name__ == "__main__":
    postprocess_HL23PV()
