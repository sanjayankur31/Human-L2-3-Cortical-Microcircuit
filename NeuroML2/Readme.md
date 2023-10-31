This folder contains versions of model components standardised in NeuroML2 format.

Learn more about NeuroML here: https://docs.neuroml.org

# Files/folders and their descriptions

- tests: tests with original Neuron cell model descriptions
- channels: all channel definitions
- synapses: all synapse definitions
- rotated_cells: folder storing all rotated cells with tonic inhibition (not tracked in Git, so you may have to re-run the network creation script to regenerate these files)
- cell_plots: plots of NeuroML cell morphologies and channel properties
- mod: original Neuron mod files and plots of their properties

- cellmorph2nml.py: exports the Neuron cell morphology to NeuroML format
- postprocess_cells.py: post process the NeuroML morphology files to add biophysics to the cells
- convert_cells_to_neuroml.sh: script to run the above two python scripts
- getinfoneuroml.py: gets information about NeuroML cell models
- create_test_network.py: creates a test network to test for validation
- create_omv_tests.py: script to create OMV test templates
- create_network.py: main script for creation of network and simulation


These files are created for different network scales:

- lems_component*.xml: files including LEMS definitions for inclusion into the LEMS simulation files
- LEMS*xml: LEMS simulation scripts
- HL23Net*.net.nml: Network description files

## NeuroML files are stored in a git-annex repository

Because of the large number of files and large sizes of the NeuroML files, these are not stored on the GitHub repository.
Instead, they are stored on [g-node](https://gin.g-node.org/ankursinha/Human-L2-3-Cortical-Microcircuit) which supports large files and git-annex.
If you want to obtain these, please follow these steps.

https://gin.g-node.org/ankursinha/Human-L2-3-Cortical-Microcircuit

### Install git-annex

It is a tool that allows for version control of large files, and is used by other tools like [datalad](https://handbook.datalad.org).
Installation instructions can be found at:

https://git-annex.branchable.com/install/


Clone the G-node repository:

```
git clone https://gin.g-node.org/ankursinha/Human-L2-3-Cortical-Microcircuit
git annex init
```
You can also use the `gin` client if you wish. Instructions are provided on g-node.

The files that are managed by git-annex will only be shown as links.
This way, they don't take up space on your local machine until you fetch them.
You can fetch individual files:

```
git annex get <path to file>
```

or get all of them:
```
git annex sync --content
```

If you don't need a file, you can "drop" it from your local hard disk:

```
git annex drop <path to file>
```
