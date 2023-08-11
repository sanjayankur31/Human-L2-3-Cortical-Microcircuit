#!/usr/bin/env python3
"""
Post process and add biophysics to cells.

We make any updates to the morphology, and add biophysics.

File: NeuroML2/postprocess_cells.py

Copyright 2022 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""


import random

import neuroml
import numpy
from neuroml.loaders import read_neuroml2_file
from neuroml.neuro_lex_ids import neuro_lex_ids
from pyneuroml.analysis import generate_current_vs_frequency_curve
from pyneuroml.lems.LEMSSimulation import LEMSSimulation
from pyneuroml.plot.PlotMorphology import plot_2D
from pyneuroml.pynml import run_lems_with_jneuroml_neuron, write_neuroml2_file

random.seed(1412)


def load_and_setup_cell(cellname: str):
    """Load a cell, and clean it to prepare it for further modifications.

    These operations are common for all cells.

    :param cellname: name of cell.
        the file containing the cell should then be <cell>.morph.cell.nml
    :returns: document with cell
    :rtype: neuroml.NeuroMLDocument

    """
    celldoc = read_neuroml2_file(f"{cellname}.morph.cell.nml")  # type: neuroml.NeuroMLDocument
    cell = celldoc.cells[0]  # type: neuroml.Cell
    celldoc.networks = []
    cell.id = cellname
    cell.notes = cell.notes.replace("NeuronTemplate_0_0", cellname)
    cell.notes += ". Reference: Yao, H. K.; Guet-McCreight, A.; Mazza, F.; Moradi Chameh, H.; Prevot, T. D.; Griffiths, J. D.; Tripathy, S. J.; Valiante, T. A.; Sibille, E. & Hay, E.  Reduced inhibition in depression impairs stimulus processing in human cortical microcircuits Cell Reports, Elsevier, 2022, 38"

    # create default groups if they don't exist
    [default_all_group, default_soma_group, default_dendrite_group, default_axon_group] = cell.setup_default_segment_groups(
        use_convention=True, default_groups=["all", "soma_group", "dendrite_group", "axon_group"]
    )

    # populate default groups
    for sg in cell.morphology.segment_groups:
        if "soma" in sg.id and sg.id != "soma_group":
            default_soma_group.add(neuroml.Include(segment_groups=sg.id))
        if "axon" in sg.id and sg.id != "axon_group":
            default_axon_group.add(neuroml.Include(segment_groups=sg.id))
        if "dend" in sg.id and sg.id != "dendrite_group":
            default_dendrite_group.add(neuroml.Include(segment_groups=sg.id))

    cell.optimise_segment_groups()

    return celldoc


def postprocess_HL23PYR():
    """Post process HL23PYR and add biophysics.

    """
    cellname = "HL23PYR"
    celldoc = load_and_setup_cell(cellname)
    cell = celldoc.cells[0]  # type: neuroml.Cell

    # apical dendrites are in groups called apic_
    # basal dendrites are in groups called dend_
    # populate the complete dendrite group, and new groups for all apical and
    # basal dendrites
    default_dendrite_group = cell.get_segment_group("dendrite_group")
    basal_group = cell.add_segment_group("basal_dendrite_group", neuro_lex_id=neuro_lex_ids["dend"], notes="Basal dendrites")
    apical_group = cell.add_segment_group("apical_dendrite_group", neuro_lex_id=neuro_lex_ids["dend"], notes="Apical dendrite_group")
    # create new global myelin group
    myelin_group = cell.add_segment_group("myelin_group", notes="Myelin group")
    for sg in cell.morphology.segment_groups:
        if "apic_" in sg.id:
            apical_group.add(neuroml.Include(segment_groups=sg.id))
        if "dend_" in sg.id:
            basal_group.add(neuroml.Include(segment_groups=sg.id))
        if "myelin_" in sg.id and sg.id != "myelin_group":
            myelin_group.add(neuroml.Include(segment_groups=sg.id))

    # myelin groups are not included in the all segment group when adding
    # biophys
    all_segment_group = cell.get_segment_group("all") # type neuroml.SegmentGroup
    all_minus_myelin = cell.add_segment_group("all_minus_myelin", notes="All without myelin")
    for inc in all_segment_group.includes:
        if "myelin" not in inc.segment_groups:
            all_minus_myelin.includes.append(inc)
    for sg in all_segment_group.members :
        if "myelin" not in sg.segment:
            all_minus_myelin.members.append(sg)

    # optimise dendrite group
    default_dendrite_group.includes = []
    default_dendrite_group.includes.append(neuroml.Include(segment_groups=apical_group.id))
    default_dendrite_group.includes.append(neuroml.Include(segment_groups=basal_group.id))

    cell.optimise_segment_groups()
    cell.reorder_segment_groups()

    # biophysics
    # include calcium dynamics component

    celldoc.add(neuroml.IncludeType(href="CaDynamics_E2_NML2.nml"), validate=False)
    # all
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="pas",
                             ion_channel="pas",
                             cond_density="0.0000954 S_per_cm2",
                             erev="-80 mV",
                             group_id="all_minus_myelin",
                             ion="non_specific",
                             ion_chan_def_file="channels/pas.channel.nml")
    cell.set_resistivity("0.1 kohm_cm", group_id="all_minus_myelin")
    cell.set_specific_capacitance("1 uF_per_cm2", group_id="all_minus_myelin")
    cell.set_init_memb_potential("-80mV")

    # myelin
    cell.set_specific_capacitance("0.02 uF_per_cm2",
                                  group_id="myelin_group")
    cell.set_resistivity("0.1 kohm_cm", group_id="myelin_group")

    # write passive cell
    write_neuroml2_file(celldoc, f"{cellname}.pas.cell.nml")

    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Ih",
                             ion_channel="Ih",
                             cond_density="0.000148 S_per_cm2",
                             erev="-45 mV",
                             group_id="all_minus_myelin",
                             ion="hcn",
                             ion_chan_def_file="channels/Ih.channel.nml")

    # somatic
    soma_group = cell.get_segment_group("soma_group")
    sgid = soma_group.id
    print(f"Adding channels to {sgid}")
    # K
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="SK_somatic",
                             ion_channel="SK",
                             cond_density="0.000853 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/SK.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="K_T_somatic",
                             ion_channel="K_T",
                             cond_density="0.0605 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/K_T.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="K_P_somatic",
                             ion_channel="K_P",
                             cond_density="0.000208 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/K_P.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Kv3_1_somatic",
                             ion_channel="Kv3_1",
                             cond_density="0.0424 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/Kv3_1.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Ih_somatic",
                             ion_channel="Ih",
                             cond_density="0.000148 S_per_cm2",
                             erev="-45 mV",
                             group_id=sgid,
                             ion="hcn",
                             ion_chan_def_file="channels/Ih.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Im_somatic",
                             ion_channel="Im",
                             cond_density="0.000306 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/Im.channel.nml")

    # Na
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="NaTg_somatic",
                             ion_channel="NaTg_PYR_somatic",
                             cond_density="0.272 S_per_cm2",
                             erev="50 mV",
                             group_id=sgid,
                             ion="na",
                             ion_chan_def_file="channels/NaTg/NaTg.channel.nml")
    # Ca
    # external concentration is set to defaults that NEURON
    # starts with
    # internal concentration is set to minCai in mod file
    cell.add_intracellular_property("Species", validate=False,
                                    id="ca",
                                    concentration_model="CaDynamics_E2_NML2_PYR_somatic",
                                    ion="ca",
                                    initial_concentration="1e-4 mM",
                                    initial_ext_concentration="2.0E-6 mol_per_cm3",
                                    segment_groups=sgid)
    # https://www.neuron.yale.edu/neuron/static/new_doc/modelspec/programmatic/ions.html
    cell.add_channel_density_v(
        "ChannelDensityNernst",
        nml_cell_doc=celldoc,
        id="Ca_HVA_somatic",
        ion_channel="Ca_HVA",
        cond_density="0.00155 S_per_cm2",
        segment_groups=sgid,
        ion="ca",
        ion_chan_def_file="channels/Ca_HVA.channel.nml")
    cell.add_channel_density_v(
        "ChannelDensityNernst",
        nml_cell_doc=celldoc,
        id="Ca_LVA_somatic",
        ion_channel="Ca_LVA",
        cond_density="0.00296 S_per_cm2",
        segment_groups=sgid,
        ion="ca",
        ion_chan_def_file="channels/Ca_LVA.channel.nml")

    # Apical
    sg = cell.get_segment_group("apical_dendrite_group")
    cell.set_specific_capacitance("2 uF_per_cm2",
                                  group_id=sg.id)
    # Add parameter that we use to distribute Ih
    sg.add(
        "InhomogeneousParameter",
        id="PathLengthOverApicDends",
        variable="p",
        metric="Path Length from root",
        proximal=sg.component_factory(
            "ProximalDetails",
            translation_start="0")
    )
    # distribute Ih
    cdnonuniform = cell.add_channel_density_v(
        "ChannelDensityNonUniform",
        nml_cell_doc=celldoc,
        id="Ih_apical",
        ion_channel="Ih",
        ion="hcn",
        erev="-45 mV",
        validate=False
    )
    varparam = cdnonuniform.add(
        "VariableParameter",
        parameter="condDensity",
        segment_groups=sg.id,
        validate=False
    )
    # TODO: clarify unit conversions
    # 606.3464 is the value of getLongestBranch("apic")
    # run test_HL23PYR.hoc, and then run:
    # > access HL23PYR.apic
    # > HL23PYR.soma distance()
    # > HL23PYR.getLongestBranch("apic")
    inhomogeneous_value = varparam.add(
        "InhomogeneousValue",
        inhomogeneous_parameters="PathLengthOverApicDends",
        value="1E9 * (0.148 * 1E-8) * ((2.087 * exp( 3.6161 * (p/606.3464))) - 0.8696)"
    )

    # Basal
    basal_group = cell.get_segment_group("basal_dendrite_group")
    sgid = basal_group.id
    cell.set_specific_capacitance("2 uF_per_cm2",
                                  group_id=sgid)
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Ih_basal",
                             ion_channel="Ih",
                             cond_density="0.000000709 S_per_cm2",
                             erev="-45 mV",
                             group_id=sgid,
                             ion="hcn",
                             ion_chan_def_file="channels/Ih.channel.nml")

    # Axonal
    axon_group = cell.get_segment_group("axon_group")
    sgid = axon_group.id
    print(f"Adding channels to {sgid}")
    # K
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="SK_axonal",
                             ion_channel="SK",
                             cond_density="0.0145 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/SK.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="K_T_axonal",
                             ion_channel="K_T",
                             cond_density="0.0424 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/K_T.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="K_P_axonal",
                             ion_channel="K_P",
                             cond_density="0.338 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/K_P.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Im_axonal",
                             ion_channel="Im",
                             cond_density="0.000 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/Im.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Kv3_1_axonal",
                             ion_channel="Kv3_1",
                             cond_density="0.941 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/Kv3_1.channel.nml")
    # Na
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="NaTg_axonal",
                             ion_channel="NaTg_PYR_axonal",
                             cond_density="1.38 S_per_cm2",
                             erev="50 mV",
                             group_id=sgid,
                             ion="na",
                             ion_chan_def_file="channels/NaTg/NaTg.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Nap_axonal",
                             ion_channel="Nap",
                             cond_density="0.00842 S_per_cm2",
                             erev="50 mV",
                             group_id=sgid,
                             ion="na",
                             ion_chan_def_file="channels/Nap.channel.nml")

    # Ca
    cell.add_channel_density_v(
        "ChannelDensityNernst",
        nml_cell_doc=celldoc,
        id="Ca_HVA_axonal",
        ion_channel="Ca_HVA",
        cond_density="0.0003060000 S_per_cm2",
        segment_groups=sgid,
        ion="ca",
        ion_chan_def_file="channels/Ca_HVA.channel.nml")
    cell.add_channel_density_v(
        "ChannelDensityNernst",
        nml_cell_doc=celldoc,
        id="Ca_LVA_axonal",
        ion_channel="Ca_LVA",
        cond_density="0.0439000000 S_per_cm2",
        segment_groups=sgid,
        ion="ca",
        ion_chan_def_file="channels/Ca_LVA.channel.nml")
    # external concentration is set to defaults that NEURON
    # starts with
    # internal concentration is set to minCai in mod file
    cell.add_intracellular_property("Species", validate=False,
                                    id="ca",
                                    concentration_model="CaDynamics_E2_NML2_PYR_axonal",
                                    ion="ca",
                                    initial_concentration="1e-4 mM",
                                    initial_ext_concentration="2.0E-6 mol_per_cm3",
                                    segment_groups=sgid)

    # L1 validation
    # cell.validate(recursive=True)
    # cell.summary(morph=False, biophys=True)
    # use pynml writer to also run L2 validation
    write_neuroml2_file(celldoc, f"{cellname}.cell.nml")


def postprocess_HL23PV():
    """Post process HL23PV and add biophysics.

    Each cell needs its biophysics to be added, so we do each cell separately.
    """
    cellname = "HL23PV"
    celldoc = load_and_setup_cell(cellname)
    cell = celldoc.cells[0]  # type: neuroml.Cell

    print("Segment groups:")
    for sg in cell.morphology.segment_groups:
        print(f"** {sg.id} **")
        print(cell.get_all_segments_in_group(sg.id))
        print()

    # biophysics
    # include calcium dynamics component
    celldoc.add(neuroml.IncludeType(href="CaDynamics_E2_NML2.nml"), validate=False)

    # all
    # Note: ion for passive channel must be "non_specific" for NEURON
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="pas",
                             ion_channel="pas",
                             cond_density="0.00011830111773572024 S_per_cm2",
                             erev="-83.92924122901199 mV",
                             group_id="all",
                             ion="non_specific",
                             ion_chan_def_file="channels/pas.channel.nml")
    cell.set_resistivity("0.1 kohm_cm", group_id="all")
    cell.set_specific_capacitance("2 uF_per_cm2", group_id="all")
    cell.set_init_memb_potential("-80mV")

    # write passive cell before adding active properties
    write_neuroml2_file(celldoc, f"{cellname}.pas.cell.nml")

    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Ih",
                             ion_channel="Ih",
                             cond_density="2.7671764064314368e-05 S_per_cm2",
                             erev="-45 mV",
                             group_id="all",
                             ion="hcn",
                             ion_chan_def_file="channels/Ih.channel.nml")

    # somatic
    soma_group = cell.get_segment_group("soma_group")
    sgid = soma_group.id
    print(f"Adding channels to {sgid}")
    # Na
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="NaTg_somatic",
                             ion_channel="NaTg_PV",
                             cond_density="0.49958525078702043 S_per_cm2",
                             erev="50 mV",
                             group_id=sgid,
                             ion="na",
                             ion_chan_def_file="channels/NaTg/NaTg.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Nap_somatic",
                             ion_channel="Nap",
                             cond_density="0.008795461417521086 S_per_cm2",
                             erev="50 mV",
                             group_id=sgid,
                             ion="na",
                             ion_chan_def_file="channels/Nap.channel.nml")

    # K
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="K_P_somatic",
                             ion_channel="K_P",
                             cond_density="9.606092478937705e-06 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/K_P.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="K_T_somatic",
                             ion_channel="K_T",
                             cond_density="0.0011701702607527396 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/K_T.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Kv3_1_somatic",
                             ion_channel="Kv3_1",
                             cond_density="2.9921080101237565 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/Kv3_1.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Im_somatic",
                             ion_channel="Im",
                             cond_density="0.04215865946497755 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/Im.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="SK_somatic",
                             ion_channel="SK",
                             cond_density="3.7265770903193036e-06 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/SK.channel.nml")
    # Ca
    # external concentration is set to defaults that NEURON
    # starts with
    # internal concentration is set to minCai in mod file
    cell.add_intracellular_property("Species", validate=False,
                                    id="ca",
                                    concentration_model="CaDynamics_E2_NML2_PV_somatic",
                                    ion="ca",
                                    initial_concentration="1e-4 mM",
                                    initial_ext_concentration="2.0E-6 mol_per_cm3",
                                    segment_groups=sgid)
    # https://www.neuron.yale.edu/neuron/static/new_doc/modelspec/programmatic/ions.html
    cell.add_channel_density_v(
        "ChannelDensityNernst",
        nml_cell_doc=celldoc,
        id="Ca_HVA_somatic",
        ion_channel="Ca_HVA",
        cond_density="0.00017953651378188165 S_per_cm2",
        segment_groups=sgid,
        ion="ca",
        ion_chan_def_file="channels/Ca_HVA.channel.nml")
    cell.add_channel_density_v(
        "ChannelDensityNernst",
        nml_cell_doc=celldoc,
        id="Ca_LVA_somatic",
        ion_channel="Ca_LVA",
        cond_density="0.09250008555398015 S_per_cm2",
        segment_groups=sgid,
        ion="ca",
        ion_chan_def_file="channels/Ca_LVA.channel.nml")

    # axonal
    sgs = cell.get_segment_group("axon_group")
    sgid = sgs.id
    print(f"Adding channels to {sgid}")
    # Na
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="NaTg_axonal",
                             ion_channel="NaTg_PV",
                             cond_density="0.10914576408883477 S_per_cm2",
                             erev="50 mV",
                             group_id=sgid,
                             ion="na",
                             ion_chan_def_file="channels/NaTg/NaTg.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Nap_axonal",
                             ion_channel="Nap",
                             cond_density="0.001200899579358837 S_per_cm2",
                             erev="50 mV",
                             group_id=sgid,
                             ion="na",
                             ion_chan_def_file="channels/Nap.channel.nml")

    # K
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="K_P_axonal",
                             ion_channel="K_P",
                             cond_density="0.6854776593761795 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/K_P.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="K_T_axonal",
                             ion_channel="K_T",
                             cond_density="0.07603372775662909 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/K_T.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Kv3_1_axonal",
                             ion_channel="Kv3_1",
                             cond_density="2.988867483754507 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/Kv3_1.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Im_axonal",
                             ion_channel="Im",
                             cond_density="0.029587905136596156 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/Im.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="SK_axonal",
                             ion_channel="SK",
                             cond_density="0.5121938998281017 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/SK.channel.nml")
    # Ca
    cell.add_channel_density_v(
        "ChannelDensityNernst",
        nml_cell_doc=celldoc,
        id="Ca_HVA_axonal",
        ion_channel="Ca_HVA",
        cond_density="0.002961469262723619 S_per_cm2",
        segment_groups=sgid,
        ion="ca",
        ion_chan_def_file="channels/Ca_HVA.channel.nml")
    cell.add_channel_density_v(
        "ChannelDensityNernst",
        nml_cell_doc=celldoc,
        id="Ca_LVA_axonal",
        ion_channel="Ca_LVA",
        cond_density="5.9457835817342756e-05 S_per_cm2",
        segment_groups=sgid,
        ion="ca",
        ion_chan_def_file="channels/Ca_LVA.channel.nml")
    # external concentration is set to defaults that NEURON
    # starts with
    # internal concentration is set to minCai in mod file
    cell.add_intracellular_property("Species", validate=False,
                                    id="ca",
                                    concentration_model="CaDynamics_E2_NML2_PV_axonal",
                                    ion="ca",
                                    initial_concentration="1e-4 mM",
                                    initial_ext_concentration="2.0E-6 mol_per_cm3",
                                    segment_groups=sgid)
    # L1 validation
    cell.validate(recursive=True)
    cell.summary(morph=True, biophys=True)
    # use pynml writer to also run L2 validation
    write_neuroml2_file(celldoc, f"{cellname}.cell.nml")


def postprocess_HL23SST():
    """Post process HL23SST and add biophysics.

    Each cell needs its biophysics to be added, so we do each cell separately.
    """
    cellname = "HL23SST"
    celldoc = load_and_setup_cell(cellname)
    cell = celldoc.cells[0]  # type: neuroml.Cell

    # basal dendrites are in groups called dend_
    # populate the complete dendrite group, and new groups for all apical and
    # basal dendrites
    default_dendrite_group = cell.get_segment_group("dendrite_group")
    basal_group = cell.add_segment_group("basal_dendrite_group", neuro_lex_id=neuro_lex_ids["dend"], notes="Basal dendrites")
    # create new global myelin group
    myelin_group = cell.add_segment_group("myelin_group", notes="Myelin group")
    for sg in cell.morphology.segment_groups:
        if "dend_" in sg.id:
            basal_group.add(neuroml.Include(segment_groups=sg.id))
        if "myelin_" in sg.id and sg.id != "myelin_group":
            myelin_group.add(neuroml.Include(segment_groups=sg.id))

    # myelin groups are not included in the all segment group when adding
    # biophys
    all_segment_group = cell.get_segment_group("all") # type neuroml.SegmentGroup
    all_minus_myelin = cell.add_segment_group("all_minus_myelin", notes="All without myelin")
    for inc in all_segment_group.includes:
        if "myelin" not in inc.segment_groups:
            all_minus_myelin.includes.append(inc)
    for sg in all_segment_group.members :
        if "myelin" not in sg.segment:
            all_minus_myelin.members.append(sg)

    # optimise dendrite group
    default_dendrite_group.includes = []
    default_dendrite_group.includes.append(neuroml.Include(segment_groups=basal_group.id))

    cell.optimise_segment_groups()
    cell.reorder_segment_groups()

    print("Segment groups:")
    for sg in cell.morphology.segment_groups:
        print(f"** {sg.id} **")
        print(cell.get_all_segments_in_group(sg.id))
        print()

    # biophysics
    # include calcium dynamics component
    celldoc.add(neuroml.IncludeType(href="CaDynamics_E2_NML2.nml"), validate=False)

    # all
    # Note: ion for passive channel must be "non_specific" for NEURON
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="pas",
                             ion_channel="pas",
                             cond_density="0.0000232 S_per_cm2",
                             erev="-81.5 mV",
                             group_id="all_minus_myelin",
                             ion="non_specific",
                             ion_chan_def_file="channels/pas.channel.nml")
    cell.set_resistivity("0.1 kohm_cm", group_id="all_minus_myelin")
    cell.set_specific_capacitance("1 uF_per_cm2", group_id="all_minus_myelin")
    cell.set_init_memb_potential("-80mV")

    # myelin
    cell.set_specific_capacitance("0.02 uF_per_cm2",
                                  group_id="myelin_group")
    cell.set_resistivity("0.1 kohm_cm", group_id="myelin_group")
    # write passive cell
    write_neuroml2_file(celldoc, f"{cellname}.pas.cell.nml")

    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Ih",
                             ion_channel="Ih",
                             cond_density="0.00001 S_per_cm2",
                             erev="-45 mV",
                             group_id="all_minus_myelin",
                             ion="hcn",
                             ion_chan_def_file="channels/Ih.channel.nml")

    # somatic
    soma_group = cell.get_segment_group("soma_group")
    sgid = soma_group.id
    print(f"Adding channels to {sgid}")
    # Na
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="NaTg_somatic",
                             ion_channel="NaTg_SST_somatic",
                             cond_density="0.127 S_per_cm2",
                             erev="50 mV",
                             group_id=sgid,
                             ion="na",
                             ion_chan_def_file="channels/NaTg/NaTg.channel.nml")

    # hcn
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Ih_somatic",
                             ion_channel="Ih",
                             cond_density="0.0000431 S_per_cm2",
                             erev="-45 mV",
                             group_id=sgid,
                             ion="hcn",
                             ion_chan_def_file="channels/Ih.channel.nml")
    # K
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="K_P_somatic",
                             ion_channel="K_P",
                             cond_density="0.0111 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/K_P.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="K_T_somatic",
                             ion_channel="K_T",
                             cond_density="0.00000 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/K_T.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Kv3_1_somatic",
                             ion_channel="Kv3_1",
                             cond_density="0.871 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/Kv3_1.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Im_somatic",
                             ion_channel="Im",
                             cond_density="0.000158 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/Im.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="SK_somatic",
                             ion_channel="SK",
                             cond_density="0.00 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/SK.channel.nml")
    # Ca
    # external concentration is set to defaults that NEURON
    # starts with
    # internal concentration is set to minCai in mod file
    cell.add_intracellular_property("Species", validate=False,
                                    id="ca",
                                    concentration_model="CaDynamics_E2_NML2_SST_somatic",
                                    ion="ca",
                                    initial_concentration="1e-4 mM",
                                    initial_ext_concentration="2.0E-6 mol_per_cm3",
                                    segment_groups=sgid)
    # https://www.neuron.yale.edu/neuron/static/new_doc/modelspec/programmatic/ions.html
    cell.add_channel_density_v(
        "ChannelDensityNernst",
        nml_cell_doc=celldoc,
        id="Ca_HVA_somatic",
        ion_channel="Ca_HVA",
        cond_density="0.00355 S_per_cm2",
        segment_groups=sgid,
        ion="ca",
        ion_chan_def_file="channels/Ca_HVA.channel.nml")
    cell.add_channel_density_v(
        "ChannelDensityNernst",
        nml_cell_doc=celldoc,
        id="Ca_LVA_somatic",
        ion_channel="Ca_LVA",
        cond_density="0.00314 S_per_cm2",
        segment_groups=sgid,
        ion="ca",
        ion_chan_def_file="channels/Ca_LVA.channel.nml")

    # Basal
    basal_group = cell.get_segment_group("basal_dendrite_group")
    sgid = basal_group.id
    cell.set_specific_capacitance("1 uF_per_cm2",
                                  group_id=sgid)
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Ih_basal",
                             ion_channel="Ih",
                             cond_density="0.0000949 S_per_cm2",
                             erev="-45 mV",
                             group_id=sgid,
                             ion="hcn",
                             ion_chan_def_file="channels/Ih.channel.nml")

    # axonal
    sgs = cell.get_segment_group("axon_group")
    sgid = sgs.id
    print(f"Adding channels to {sgid}")
    # Na
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="NaTg_axonal",
                             ion_channel="NaTg_SST_axonal",
                             cond_density="0.343 S_per_cm2",
                             erev="50 mV",
                             group_id=sgid,
                             ion="na",
                             ion_chan_def_file="channels/NaTg/NaTg.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Nap_axonal",
                             ion_channel="Nap",
                             cond_density="0.000444 S_per_cm2",
                             erev="50 mV",
                             group_id=sgid,
                             ion="na",
                             ion_chan_def_file="channels/Nap.channel.nml")

    # K
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="K_P_axonal",
                             ion_channel="K_P",
                             cond_density="0.0295 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/K_P.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="K_T_axonal",
                             ion_channel="K_T",
                             cond_density="0.023 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/K_T.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Kv3_1_axonal",
                             ion_channel="Kv3_1",
                             cond_density="0.984 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/Kv3_1.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Im_axonal",
                             ion_channel="Im",
                             cond_density="0.000317 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/Im.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="SK_axonal",
                             ion_channel="SK",
                             cond_density="0.00113 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/SK.channel.nml")
    # Ca
    cell.add_channel_density_v(
        "ChannelDensityNernst",
        nml_cell_doc=celldoc,
        id="Ca_HVA_axonal",
        ion_channel="Ca_HVA",
        cond_density="0.00145 S_per_cm2",
        segment_groups=sgid,
        ion="ca",
        ion_chan_def_file="channels/Ca_HVA.channel.nml")
    cell.add_channel_density_v(
        "ChannelDensityNernst",
        nml_cell_doc=celldoc,
        id="Ca_LVA_axonal",
        ion_channel="Ca_LVA",
        cond_density="0.0627 S_per_cm2",
        segment_groups=sgid,
        ion="ca",
        ion_chan_def_file="channels/Ca_LVA.channel.nml")
    # external concentration is set to defaults that NEURON
    # starts with
    # internal concentration is set to minCai in mod file
    cell.add_intracellular_property("Species", validate=False,
                                    id="ca",
                                    concentration_model="CaDynamics_E2_NML2_SST_axonal",
                                    ion="ca",
                                    initial_concentration="1e-4 mM",
                                    initial_ext_concentration="2.0E-6 mol_per_cm3",
                                    segment_groups=sgid)
    # L1 validation
    cell.validate(recursive=True)
    cell.summary(morph=True, biophys=True)
    # use pynml writer to also run L2 validation
    write_neuroml2_file(celldoc, f"{cellname}.cell.nml")


def postprocess_HL23VIP():
    """Post process HL23VIP and add biophysics.

    Each cell needs its biophysics to be added, so we do each cell separately.
    """
    cellname = "HL23VIP"
    celldoc = load_and_setup_cell(cellname)
    cell = celldoc.cells[0]  # type: neuroml.Cell

    print("Segment groups:")
    for sg in cell.morphology.segment_groups:
        print(f"** {sg.id} **")
        print(cell.get_all_segments_in_group(sg.id))
        print()

    # biophysics
    # include calcium dynamics component
    celldoc.add(neuroml.IncludeType(href="CaDynamics_E2_NML2.nml"), validate=False)

    # all
    # Note: ion for passive channel must be "non_specific" for NEURON
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="pas",
                             ion_channel="pas",
                             cond_density="2.5756438955642182e-05 S_per_cm2",
                             erev="-79.74132024971513 mV",
                             group_id="all",
                             ion="non_specific",
                             ion_chan_def_file="channels/pas.channel.nml")
    cell.set_resistivity("0.1 kohm_cm", group_id="all")
    cell.set_specific_capacitance("2 uF_per_cm2", group_id="all")
    cell.set_init_memb_potential("-80mV")

    # write passive cell
    write_neuroml2_file(celldoc, f"{cellname}.pas.cell.nml")

    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Ih",
                             ion_channel="Ih",
                             cond_density="4.274951616063423e-05 S_per_cm2",
                             erev="-45 mV",
                             group_id="all",
                             ion="hcn",
                             ion_chan_def_file="channels/Ih.channel.nml")

    # somatic
    soma_group = cell.get_segment_group("soma_group")
    sgid = soma_group.id
    print(f"Adding channels to {sgid}")
    # Na
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="NaTg_somatic",
                             ion_channel="NaTg_VIP_somatic",
                             cond_density="0.11491205828369114 S_per_cm2",
                             erev="50 mV",
                             group_id=sgid,
                             ion="na",
                             ion_chan_def_file="channels/NaTg/NaTg.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Nap_somatic",
                             ion_channel="Nap",
                             cond_density="0.0001895305240694194 S_per_cm2",
                             erev="50 mV",
                             group_id=sgid,
                             ion="na",
                             ion_chan_def_file="channels/Nap.channel.nml")

    # K
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="K_P_somatic",
                             ion_channel="K_P",
                             cond_density="0.0009925418924114282 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/K_P.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="K_T_somatic",
                             ion_channel="K_T",
                             cond_density="0.009051981253674193 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/K_T.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Kv3_1_somatic",
                             ion_channel="Kv3_1",
                             cond_density="0.31215653649208114 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/Kv3_1.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Im_somatic",
                             ion_channel="Im",
                             cond_density="0.0003679378262289559 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/Im.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="SK_somatic",
                             ion_channel="SK",
                             cond_density="0.1655502166633749 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/SK.channel.nml")
    # Ca
    # external concentration is set to defaults that NEURON
    # starts with
    # internal concentration is set to minCai in mod file
    cell.add_intracellular_property("Species", validate=False,
                                    id="ca",
                                    concentration_model="CaDynamics_E2_NML2_VIP_somatic",
                                    ion="ca",
                                    initial_concentration="1e-4 mM",
                                    initial_ext_concentration="2.0E-6 mol_per_cm3",
                                    segment_groups=sgid)
    # https://www.neuron.yale.edu/neuron/static/new_doc/modelspec/programmatic/ions.html
    cell.add_channel_density_v(
        "ChannelDensityNernst",
        nml_cell_doc=celldoc,
        id="Ca_HVA_somatic",
        ion_channel="Ca_HVA",
        cond_density="4.384846294634834e-05 S_per_cm2",
        segment_groups=sgid,
        ion="ca",
        ion_chan_def_file="channels/Ca_HVA.channel.nml")
    cell.add_channel_density_v(
        "ChannelDensityNernst",
        nml_cell_doc=celldoc,
        id="Ca_LVA_somatic",
        ion_channel="Ca_LVA",
        cond_density="0.0034472458995879864 S_per_cm2",
        segment_groups=sgid,
        ion="ca",
        ion_chan_def_file="channels/Ca_LVA.channel.nml")

    # axonal
    sgs = cell.get_segment_group("axon_group")
    sgid = sgs.id
    print(f"Adding channels to {sgid}")
    # Na
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="NaTg_axonal",
                             ion_channel="NaTg_VIP_axonal",
                             cond_density="0.20112200814143477 S_per_cm2",
                             erev="50 mV",
                             group_id=sgid,
                             ion="na",
                             ion_chan_def_file="channels/NaTg/NaTg.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Nap_axonal",
                             ion_channel="Nap",
                             cond_density="0.0006248906854665301 S_per_cm2",
                             erev="50 mV",
                             group_id=sgid,
                             ion="na",
                             ion_chan_def_file="channels/Nap.channel.nml")

    # K
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="K_P_axonal",
                             ion_channel="K_P",
                             cond_density="0.26489876414660096 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/K_P.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="K_T_axonal",
                             ion_channel="K_T",
                             cond_density="0.014364427062274185 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/K_T.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Kv3_1_axonal",
                             ion_channel="Kv3_1",
                             cond_density="0.0011201608191112877 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/Kv3_1.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Im_axonal",
                             ion_channel="Im",
                             cond_density="0.00013891465461042372 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/Im.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="SK_axonal",
                             ion_channel="SK",
                             cond_density="0.7027792087501376 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/SK.channel.nml")
    # Ca
    cell.add_channel_density_v(
        "ChannelDensityNernst",
        nml_cell_doc=celldoc,
        id="Ca_HVA_axonal",
        ion_channel="Ca_HVA",
        cond_density="2.819397237794038e-05 S_per_cm2",
        segment_groups=sgid,
        ion="ca",
        ion_chan_def_file="channels/Ca_HVA.channel.nml")
    cell.add_channel_density_v(
        "ChannelDensityNernst",
        nml_cell_doc=celldoc,
        id="Ca_LVA_axonal",
        ion_channel="Ca_LVA",
        cond_density="0.010354001513952075 S_per_cm2",
        segment_groups=sgid,
        ion="ca",
        ion_chan_def_file="channels/Ca_LVA.channel.nml")
    # external concentration is set to defaults that NEURON
    # starts with
    # internal concentration is set to minCai in mod file
    cell.add_intracellular_property("Species", validate=False,
                                    id="ca",
                                    concentration_model="CaDynamics_E2_NML2_VIP_axonal",
                                    ion="ca",
                                    initial_concentration="1e-4 mM",
                                    initial_ext_concentration="2.0E-6 mol_per_cm3",
                                    segment_groups=sgid)

    # L1 validation
    cell.validate(recursive=True)
    cell.summary(morph=True, biophys=True)
    # use pynml writer to also run L2 validation
    write_neuroml2_file(celldoc, f"{cellname}.cell.nml")


def analyse_HL23PYR(hyperpolarising: bool = True, depolarising: bool = True):
    """Generate various curves for HL23PYR cells

    :returns: None

    """
    cellname = "HL23PYR"
    if hyperpolarising:
        # hyper-polarising inputs
        generate_current_vs_frequency_curve(
            nml2_file=f"{cellname}.cell.nml",
            cell_id=cellname,
            custom_amps_nA=list(numpy.arange(-0.05, -0.1, -0.01)),
            temperature="34 degC",
            pre_zero_pulse=200,
            post_zero_pulse=300,
            plot_voltage_traces=True,
            plot_iv=True,
            plot_if=False,
            simulator="jNeuroML_NEURON",
            analysis_delay=300.,
            analysis_duration=400.
        )

    if depolarising:
        # depolarising inputs
        generate_current_vs_frequency_curve(
            nml2_file=f"{cellname}.cell.nml",
            cell_id=cellname,
            plot_voltage_traces=True,
            spike_threshold_mV=-10.0,
            # custom_amps_nA=list(numpy.arange(0, 0.3, 0.05)),
            custom_amps_nA=[0.2],
            temperature="34 degC",
            pre_zero_pulse=200,
            post_zero_pulse=300,
            plot_iv=True,
            simulator="jNeuroML_NEURON",
            analysis_delay=50.,
            analysis_duration=200.
        )


def analyse_HL23PV(hyperpolarising: bool = True, depolarising: bool = True):
    """Generate various curves for HL23PV cells

    :returns: None

    """
    cellname = "HL23PV"
    if hyperpolarising:
        # hyper-polarising inputs
        generate_current_vs_frequency_curve(
            nml2_file=f"{cellname}.cell.nml",
            cell_id=cellname,
            custom_amps_nA=list(numpy.arange(-0.05, -0.1, -0.01)),
            pre_zero_pulse=200,
            post_zero_pulse=300,
            plot_voltage_traces=True,
            plot_iv=True,
            plot_if=False,
            simulator="jNeuroML_NEURON",
            analysis_delay=300.,
            analysis_duration=400.
        )

    if depolarising:
        # depolarising inputs
        generate_current_vs_frequency_curve(
            nml2_file=f"{cellname}.cell.nml",
            cell_id=cellname,
            plot_voltage_traces=True,
            spike_threshold_mV=-10.0,
            custom_amps_nA=list(numpy.arange(0, 0.3, 0.01)),
            pre_zero_pulse=200,
            post_zero_pulse=300,
            plot_iv=True,
            simulator="jNeuroML_NEURON",
            analysis_delay=300.,
            analysis_duration=400.
        )


def analyse_HL23SST(hyperpolarising: bool = True, depolarising: bool = True):
    """Generate various curves for HL23SST cells

    :returns: None

    """
    cellname = "HL23SST"
    if hyperpolarising:
        # hyper-polarising inputs
        generate_current_vs_frequency_curve(
            nml2_file=f"{cellname}.cell.nml",
            cell_id=cellname,
            custom_amps_nA=list(numpy.arange(-0.05, -0.1, -0.01)),
            pre_zero_pulse=200,
            post_zero_pulse=300,
            plot_voltage_traces=True,
            plot_iv=True,
            plot_if=False,
            simulator="jNeuroML_NEURON",
            analysis_delay=300.,
            analysis_duration=400.
        )

    if depolarising:
        # depolarising inputs
        generate_current_vs_frequency_curve(
            nml2_file=f"{cellname}.cell.nml",
            cell_id=cellname,
            plot_voltage_traces=True,
            spike_threshold_mV=-10.0,
            # custom_amps_nA=list(numpy.arange(0, 0.3, 0.01)),
            custom_amps_nA=[0.2],
            pre_zero_pulse=200,
            post_zero_pulse=300,
            plot_iv=True,
            simulator="jNeuroML_NEURON",
            analysis_delay=300.,
            analysis_duration=400.
        )


def analyse_HL23VIP(hyperpolarising: bool = True, depolarising: bool = True):
    """Generate various curves for HL23VIP cells

    :returns: None

    """
    cellname = "HL23VIP"
    if hyperpolarising:
        # hyper-polarising inputs
        generate_current_vs_frequency_curve(
            nml2_file=f"{cellname}.cell.nml",
            cell_id=cellname,
            custom_amps_nA=list(numpy.arange(-0.05, -0.1, -0.01)),
            pre_zero_pulse=200,
            post_zero_pulse=300,
            plot_voltage_traces=True,
            plot_iv=True,
            plot_if=False,
            simulator="jNeuroML_NEURON",
            analysis_delay=300.,
            analysis_duration=400.
        )

    if depolarising:
        # depolarising inputs
        generate_current_vs_frequency_curve(
            nml2_file=f"{cellname}.cell.nml",
            cell_id=cellname,
            plot_voltage_traces=True,
            spike_threshold_mV=-10.0,
            # custom_amps_nA=list(numpy.arange(0, 0.3, 0.01)),
            custom_amps_nA=[0.2],
            pre_zero_pulse=200,
            post_zero_pulse=300,
            plot_iv=True,
            simulator="jNeuroML_NEURON",
            analysis_delay=300.,
            analysis_duration=400.
        )


def simulate_test_network(cells: list = []):
    """Create and simulate a test network

    WIP: does not simulate it yet.

    :param cells: list of cell names to simulate
    :type cells: list of strings
    :returns: None

    """
    net_doc = neuroml.NeuroMLDocument.component_factory("NeuroMLDocument",
                                                        id="HL23")

    network = net_doc.add("Network", id="HL23Net",
                          type="networkWithTemperature", temperature="34 degC", validate=False)  # type: neuroml.Network

    for cell in cells:
        # a document for the cell with input for OMV test
        cell_net_doc = neuroml.NeuroMLDocument.component_factory("NeuroMLDocument", id=f"{cell}")
        cell_network = cell_net_doc.add("Network", id=f"{cell}Net",
                                        type="networkWithTemperature", temperature="34 degC", validate=False)  # type: neuroml.Network
        cell_net_doc.add("IncludeType", href=f"{cell}.cell.nml")
        cell_pop = cell_network.add("Population", id=f"{cell}_pop",
                                    type="populationList", component=f"{cell}",
                                    validate=False)  # type: neuroml.Population
        cell_pop.add("Property", tag="color", value=f"{random.random()} {random.random()} {random.random()}")
        cell_pop.add("Property", tag="region", value="L23")
        cell_pop.add(
            "Instance", id="0",
            location=cell_pop.component_factory(
                "Location", x="0", y="0", z="0",
            ))
        # to match the test_*hoc NEURON files
        cell_net_doc.add("PulseGenerator", id=f"pg_{cell}", notes="Simple pulse generator", delay="50ms", duration="200ms", amplitude="0.2nA")
        cell_inputs = cell_network.add("InputList", id=f"stim_iclamp_{cell}",
                                       populations=f"{cell}_pop", component=f"pg_{cell}",
                                       validate=False)
        cell_inputs.add("Input", id=0, target=f"../{cell}_pop/0",
                        destination="synapses")
        cell_net_doc.validate(True)
        write_neuroml2_file(cell_net_doc, f"{cell}.net.nml")
        cell_net_sim = LEMSSimulation(f"{cell}_sim", duration=300., dt=0.01,
                                      target=cell_network.id)
        cell_net_sim.include_neuroml2_file(f"{cell}.net.nml")
        cell_net_sim.create_output_file(id=f"{cell}_net_output",
                                        file_name=f"{cell}_net.dat")
        cell_net_sim.add_column_to_output_file(output_file_id=f"{cell}_net_output",
                                               column_id=f"{cell}_pop_0_v",
                                               quantity=f"{cell}_pop/0/{cell}/0/v")
        # Save LEMS simulation to file
        cell_sim_file = cell_net_sim.save_to_file()

        # Run the simulation using the NEURON simulator
        run_lems_with_jneuroml_neuron(cell_sim_file, max_memory="8G", nogui=True,
                                      plot=False, skip_run=False)

        # network with all 4 cells
        net_doc.add("IncludeType", href=f"{cell}.cell.nml")
        pop = network.add("Population", id=f"{cell}_pop",
                          type="populationList", component=f"{cell}",
                          validate=False)  # type: neuroml.Population
        pop.add("Property", tag="color", value=f"{random.random()} {random.random()} {random.random()}")
        pop.add("Property", tag="region", value="L23")
        pop.add(
            "Instance", id="0",
            location=pop.component_factory(
                "Location",
                x=f"{random.randint(0, 200)}",
                y=f"{random.randint(0, 200)}",
                z=f"{random.randint(0, 200)}")
        )
        net_doc.add("PulseGenerator", id=f"pg_{cell}", notes="Simple pulse generator", delay="50ms", duration="200ms", amplitude="0.2nA")
        network.add("ExplicitInput", target="pop0[0]", input=f"pg_{cell}")

    net_doc.validate(True)
    write_neuroml2_file(net_doc, "HL23.net.nml")
    # plot_2D("HL23.net.nml", plane2d="xy")


if __name__ == "__main__":
    cellnames = ["HL23PV", "HL23PYR", "HL23SST", "HL23VIP"]
    
    import sys 
    if '-postprocall' in sys.argv:
        postprocess_HL23VIP()
        postprocess_HL23PV()
        postprocess_HL23PYR()
        postprocess_HL23SST()


    """
    postprocess_HL23VIP()
    analyse_HL23VIP(True, True)
    simulate_test_network(cellnames)
    postprocess_HL23PV()
    analyse_HL23PV(True, True)
    postprocess_HL23PYR()
    analyse_HL23PYR(False, True)
    postprocess_HL23SST()
    analyse_HL23SST(True, True)
    """
