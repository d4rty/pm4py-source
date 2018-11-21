import os
from pm4py.util.create_artificial_event_log import create_artificial_event_log
from pm4py.objects.log.importer.xes import factory as xes_importer
from pm4py.algo.conformance.alignments.most_probable_alignments.utils import *
from pm4py.objects.petri.petrinet import PetriNet
from pm4py.objects.petri import utils as petri_net_utils
from pm4py.visualization.petrinet import factory as petri_net_visualization_factory
from pm4py.objects.petri.petrinet import Marking

event_log_path = os.path.join('C:/', 'Users', 'Daniel', 'Desktop', 'master_thesis_code', 'pm4py-source-forked', 'tests',
                              'input_data')
event_log_name = "trace_from_paper_automatic_root_cause_identification"

traces = [
    {"frequency": 40, "events": ["A", "B", "C"]},
    {"frequency": 10, "events": ["A", "B", "D"]},
    {"frequency": 80, "events": ["A", "B", "B", "C"]},
    {"frequency": 20, "events": ["A", "B", "B", "D"]},
    {"frequency": 40, "events": ["A", "B", "B", "B", "C"]},
    {"frequency": 10, "events": ["A", "B", "B", "B", "D"]},
    {"frequency": 1, "events": ["A", "B", "B", "B", "B", "B", "B", "B", "C"]},
    {"frequency": 1, "events": ["A", "B", "B", "B", "B"]},
    {"frequency": 1, "events": ["A", "C", "B"]},
]

create_artificial_event_log(event_log_path, event_log_name, traces, force=True)
event_log = xes_importer.import_log(os.path.join(event_log_path, event_log_name) + '.xes')
prior = {'A': 1, 'B': 1, 'C': 1, 'D': 1, '*': 1}

# create petri net #####################################################################################################
petri_net = PetriNet("running_example")
places = {}
for i in range(1, 5):
    places['p_%i' % i] = PetriNet.Place('p_%i' % i)
    petri_net.places.add(places['p_%i' % i])

transitions = {
    # PetriNet.Transition(<unique name in petri net>, <label>)
    't_A': PetriNet.Transition('t_A', 'A'),
    't_B': PetriNet.Transition('t_B', 'B'),
    't_C': PetriNet.Transition('t_C', 'C'),
    't_D': PetriNet.Transition('t_D', 'D'),
    't_1': PetriNet.Transition('t_1', None),  # silent transition, since label=None
}

for transition in transitions:
    # print(transition)
    petri_net.transitions.add(transitions[transition])

petri_net_utils.add_arc_from_to(places['p_1'], transitions['t_A'], petri_net)
petri_net_utils.add_arc_from_to(transitions['t_A'], places['p_2'], petri_net)
petri_net_utils.add_arc_from_to(places['p_2'], transitions['t_B'], petri_net)
petri_net_utils.add_arc_from_to(transitions['t_B'], places['p_3'], petri_net)
petri_net_utils.add_arc_from_to(places['p_3'], transitions['t_C'], petri_net)
petri_net_utils.add_arc_from_to(places['p_3'], transitions['t_D'], petri_net)
petri_net_utils.add_arc_from_to(places['p_3'], transitions['t_1'], petri_net)
petri_net_utils.add_arc_from_to(transitions['t_1'], places['p_2'], petri_net)
petri_net_utils.add_arc_from_to(transitions['t_C'], places['p_4'], petri_net)
petri_net_utils.add_arc_from_to(transitions['t_D'], places['p_4'], petri_net)

initial_marking = Marking()
initial_marking[places['p_1']] = 1
final_marking = Marking()
final_marking[places['p_4']] = 1

# petri net visualization
parameters = {}
parameters["format"] = "svg"
# parameters['debug'] = True
# parameters = None
# gviz = petri_net_visualization_factory.apply(petri_net, initial_marking, final_marking, parameters=parameters)
# petri_net_visualization_factory.view(gviz)

# end create petri net #################################################################################################


log_move_probabilities = calculate_log_move_probability(event_log, prior)
model_move_probability = calculate_model_move_probability(event_log, petri_net, initial_marking, final_marking)
