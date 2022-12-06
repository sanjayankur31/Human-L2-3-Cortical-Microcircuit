from neuromllite import *
from neuromllite.NetworkGenerator import *
from neuromllite.utils import create_new_model
import sys


colors = {'HL23PV':'0 0.8 0', 'HL23PYR':'0.8 0 0'}

def generate(cell, duration=300, config='IClamp',parameters = None):

    reference = "%s_%s"%(config, cell)

    cell_id = '%s'%cell
    cell_nmll = Cell(id=cell_id, neuroml2_source_file='%s.cell.nml'%(cell))

    ################################################################################
    ###   Add some inputs

    if 'IClamp' in config:

        if not parameters:
            parameters = {}
            parameters['stim_amp'] = '200pA'

        input_source = InputSource(id='iclamp_0',
                                   neuroml2_input='PulseGenerator',
                                   parameters={'amplitude':'stim_amp', 'delay':'50ms', 'duration':'200ms'})


    else:

        if not parameters:
            parameters = {}
            parameters['average_rate'] = '100 Hz'
            parameters['number_per_cell'] = '10'

        input_source = InputSource(id='pfs0',
                                   neuroml2_input='PoissonFiringSynapse',
                                   parameters={'average_rate':'average_rate',
                                               'synapse':syn_exc.id,
                                               'spike_target':"./%s"%syn_exc.id})

    sim, net = create_new_model(reference,
                     duration,
                     dt=0.01, # ms
                     temperature=34.0, # degC
                     parameters = parameters,
                     cell_for_default_population=cell_nmll,
                     color_for_default_population=colors[cell],
                     input_for_default_population=input_source)

    return sim, net



if __name__ == "__main__":

    if '-all' in sys.argv:
        for cell in colors:
            sim, net = generate(cell, 300, config="IClamp",parameters={'stim_amp':'200pA'})
            check_to_generate_or_run(sys.argv, sim)


    else:

        #sim, net = generate('cADpyr229_L23_PC_c292d67a2e_0_0', 3000, config="IClamp")
        sim, net = generate('HL23PV', 300, config="IClamp",parameters={'stim_amp':'200pA'})

        check_to_generate_or_run(sys.argv, sim)
