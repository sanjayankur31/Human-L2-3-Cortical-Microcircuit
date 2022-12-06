from neuromllite import *
from neuromllite.NetworkGenerator import *
from neuromllite.utils import create_new_model
import sys


colors = {'HL23PV':'0 0.9 0', 'HL23PYR':'0.9 0 0'}

def generate(cell_numbers, duration=300, config='IClamp', parameters = None):

    reference = "%s"%(config)
    cells_nmll = {}

    for cell in cell_numbers:

        cell_id = '%s'%cell
        reference+='_%s'%cell_id
        cells_nmll[cell_id] = Cell(id=cell_id, neuroml2_source_file='%s.cell.nml'%(cell))

    ################################################################################

    if 'IClamp' in config:

        if not parameters:
            parameters = {}
            parameters['stim_amp'] = '200pA'

        input_source = InputSource(id='iclamp_0',
                                   neuroml2_input='PulseGenerator',
                                   parameters={'amplitude':'stim_amp',
                                               'delay':'50ms',
                                               'duration':'200ms'})
    else:

        if not parameters:
            parameters = {}
            parameters['average_rate'] = '100 Hz'
            parameters['number_per_cell'] = '10'

        net = Network(id=config)
        net.notes = "Network example"

        net.seed = 1234
        net.temperature = 34.0

        net.parameters = parameters

        ampa = Synapse(id="AMPA_syn", neuroml2_source_file="synapses/AMPA_syn.synapse.nml")
        net.synapses.append(ampa)

        input_source = InputSource(id='pfs0',
                                   neuroml2_input='PoissonFiringSynapse',
                                   parameters={'average_rate':'average_rate',
                                               'synapse':ampa.id,
                                               'spike_target':"./%s"%ampa.id})
        net.input_sources.append(input_source)

    if 'IClamp' in config:

        cell_nmll = list(cells_nmll.values())[0]
        sim, net = create_new_model(reference,
                         duration,
                         dt=0.01, # ms
                         temperature=34.0, # degC
                         parameters = parameters,
                         cell_for_default_population=cell_nmll,
                         color_for_default_population=colors[cell_nmll.id],
                         input_for_default_population=input_source)

    else:

        for cell in cell_numbers:
            cell_nmll = cells_nmll[cell]
            net.cells.append(cell_nmll)

            size = cell_numbers[cell]

            net.parameters['num_%s'%cell] = size

            pE = Population(
                id="Pop_%s"%cell,
                size='num_%s'%cell,
                component=cell,
                properties={"color": colors[cell]},
            )
            net.populations.append(pE)


            net.inputs.append(
                Input(
                    id="stim_%s"%cell,
                    input_source=input_source.id,
                    population=pE.id,
                    percentage=100,
                    weight=0.1,
                )
            )


        net_file = net.to_json_file("%s.json" % net.id)

        sim = Simulation(
            id="Sim_%s"%net.id,
            network=net_file,
            duration="300",
            seed="1111",
            dt="0.01",
            record_traces={"all": "*"},
        )

        sim.to_json_file()


    return sim, net



if __name__ == "__main__":

    if '-all' in sys.argv or '-net' in sys.argv:

        if '-all' in sys.argv:
            for cell in colors:
                sim, net = generate({cell:1}, 300, config="IClamp",parameters={'stim_amp':'200pA'})
                check_to_generate_or_run(sys.argv, sim)

        sim, net = generate({'HL23PV':1, 'HL23PYR':1}, 300, config="TestNetwork", parameters={'average_rate':'100 Hz'})
        check_to_generate_or_run(sys.argv, sim)

    else:

        #sim, net = generate('cADpyr229_L23_PC_c292d67a2e_0_0', 3000, config="IClamp")
        sim, net = generate({'HL23PV':1}, 300, config="IClamp",parameters={'stim_amp':'200pA'})

        check_to_generate_or_run(sys.argv, sim)
