import os

from pm4py.objects import petri
from pm4py.objects.petri.petrinet import PetriNet
from pm4py.objects.petri import utils
from pm4py.visualization.petrinet import factory as pn_vis_factory
from pm4py.objects.petri.petrinet import Marking
from pm4py.objects.log.importer.xes import factory as xes_import_factory
from pm4py.algo.conformance import alignments
from pm4py.algo.conformance.alignments.versions.state_equation_a_star import PARAM_MODEL_COST_FUNCTION
from pm4py.algo.conformance.alignments.versions.state_equation_a_star import PARAM_SYNC_COST_FUNCTION
from pm4py.algo.conformance.alignments.versions.state_equation_a_star import PARAM_TRACE_COST_FUNCTION

net = PetriNet("running_example_of_paper")

places = {}
for i in range(1, 5):
    places['p_%i' % i] = PetriNet.Place('p_%i' % i)

for place in places:
    net.places.add(places[place])

transitions = {
    'a': PetriNet.Transition('a', 'a'),
    'b': PetriNet.Transition('b', 'b'),
    'c': PetriNet.Transition('c', 'c'),
    'd': PetriNet.Transition('d', 'd'),
    't_1': PetriNet.Transition('t_1', None),  # silent transition
}

for transition in transitions:
    # print(transition)
    net.transitions.add(transitions[transition])

utils.add_arc_from_to(places['p_1'], transitions['a'], net)
utils.add_arc_from_to(transitions['a'], places['p_2'], net)
utils.add_arc_from_to(transitions['b'], places['p_3'], net)
utils.add_arc_from_to(transitions['t_1'], places['p_2'], net)
utils.add_arc_from_to(places['p_2'], transitions['b'], net)
utils.add_arc_from_to(places['p_3'], transitions['t_1'], net)
utils.add_arc_from_to(places['p_3'], transitions['c'], net)
utils.add_arc_from_to(places['p_3'], transitions['d'], net)
utils.add_arc_from_to(transitions['c'], places['p_4'], net)
utils.add_arc_from_to(transitions['d'], places['p_4'], net)

initial_marking = Marking()
initial_marking[places['p_1']] = 1
final_marking = Marking()
final_marking[places['p_4']] = 1

parameters = {}
parameters["format"] = "svg"
# parameters['debug'] = True

# parameters = None
gviz = pn_vis_factory.apply(net, initial_marking, final_marking, parameters=parameters)
# pn_vis_factory.view(gviz)

# LOG
log = xes_import_factory.apply(os.path.join("..", "tests", "input_data", "paper_most_probable_alignments.xes"))

# Alignment computation
model_cost_function = {}
sync_cost_function = {}
for t in net.transitions:
    if t.label is not None:
        model_cost_function[t] = 1
        sync_cost_function[t] = 0
    else:
        model_cost_function[t] = 1

number_synchronous_moves_per_transition = {}
number_model_moves_per_transition = {}

# assigns each trace (id) a frequency
frequency_trace = {'1': 40, '2': 10, '3': 80, '4': 20, '5': 40, '6': 10, '7': 1, '8': 1, '9': 1}

for trace in log:
    print(trace)

    # giving each event in the trace the same cost
    trace_costs = list(map(lambda event: 1, trace))

    params = {PARAM_MODEL_COST_FUNCTION: model_cost_function,
              PARAM_TRACE_COST_FUNCTION: trace_costs,
              PARAM_SYNC_COST_FUNCTION: sync_cost_function}

    alignment = alignments.factory.apply(trace, net, initial_marking, final_marking, parameters=params)

    # for simplicity, in the log file each trace occurs only once
    trace_id = trace.__getattribute__('attributes')['concept:name']
    trace_frequency = frequency_trace[trace_id]

    for move in alignment['alignment']:
        if move[1] == '>>':
            # according to paper, log moves are not considered
            continue
        # either a synchronous or model move
        if move[0] == '>>':
            # model move
            if move[1] in number_model_moves_per_transition:
                number_model_moves_per_transition[move[1]] += 1 * trace_frequency
            else:
                number_model_moves_per_transition[move[1]] = 1 * trace_frequency
        if move[0] == move[1]:
            # synchronous move
            if move[0] in number_synchronous_moves_per_transition:
                number_synchronous_moves_per_transition[move[0]] += 1 * trace_frequency
            else:
                number_synchronous_moves_per_transition[move[0]] = 1 * trace_frequency

    print("alignment length: " + str(len(alignment['alignment'])))
    print(alignment)
    print()
    print()

print("number_synchronous_moves_per_transition")
print(number_synchronous_moves_per_transition)
print("number_model_moves_per_transition (moves in model only)")
print(number_model_moves_per_transition)
print('---------------------------------------------------------')

for key, value in number_synchronous_moves_per_transition.items():
    if key in number_model_moves_per_transition:
        number_synchronous_moves_per_transition[key] += number_model_moves_per_transition[key]
print(number_synchronous_moves_per_transition)
