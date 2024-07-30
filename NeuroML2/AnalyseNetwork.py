import neuromllite
from pyneuroml import pynml
import os

from neuromllite.MatrixHandler import MatrixHandler
from neuromllite.GraphVizHandler import GraphVizHandler

f = 'HL23Net_1.0.net.nml'
#f = 'HL23Net_0.1.net.nml'

level = 1

print("Converting %s to matrix form, level %i" % (f, level))

from neuroml.hdf5.NeuroMLXMLParser import NeuroMLXMLParser

handler = MatrixHandler(level=level, 
                        nl_network=None,
                        show_already=False,
                        save_figs_to_dir='.')

currParser = NeuroMLXMLParser(handler)

currParser.parse(f)

handler.finalise_document()

level_gv = 6
engine = 'dot'

def generate_graph(level_gv, engine):
    handler_gv = GraphVizHandler(level_gv, 
                                engine=engine, 
                                nl_network=None,         
                                include_ext_inputs=False,
                                include_input_pops=False,
                                view_on_render=False)
            
    currParser_gv = NeuroMLXMLParser(handler_gv)

    currParser_gv.parse(f)

    handler_gv.finalise_document()

    os.rename('HL23Network.gv.png','%s.%s.%s.png'%(f,engine,level_gv))

generate_graph(level_gv, engine)

level_gv = 2
engine = 'circo'

generate_graph(level_gv, engine)


info = '## Analysis of NeuroML network: %s\n\n'%(f)

nml2_doc = pynml.read_neuroml2_file(f)

net_info = nml2_doc.summary(show_includes=False)
info +='```\n%s\n```\n'%(net_info)

info +='![fig](%s.circo.2.png)\n'%(f)
info +='![fig](%s.dot.6.png)\n'%(f)

for w in handler.weight_arrays_to_show:

    info +='### %s\n'%w
    info +='![fig](%s)\n'%(handler.weight_array_figures[w])
    info +='```\n%s\n```\n'%(handler.weight_arrays_to_show[w])
    #weight_array_figures


with open('Analysis_%s.md'%f,'w') as rm_file:
    rm_file.write(info)

