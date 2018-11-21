import os
from pm4py.objects.petri.petrinet import PetriNet
from pm4py.objects.petri import utils
from pm4py.visualization.petrinet import factory as pn_vis_factory
from pm4py.objects.petri.petrinet import Marking
from pm4py.objects.log.importer.xes import factory as xes_import_factory
from pm4py.algo.conformance import alignments
from pm4py.algo.conformance.alignments.versions.state_equation_a_star import PARAM_MODEL_COST_FUNCTION
from pm4py.algo.conformance.alignments.versions.state_equation_a_star import PARAM_SYNC_COST_FUNCTION
from pm4py.algo.conformance.alignments.versions.state_equation_a_star import PARAM_TRACE_COST_FUNCTION


net = PetriNet("petri_net_m_1_paper_replaying_history")

places = {}
for i in range(1, 8):
    places['p_%i' % i] = PetriNet.Place('p_%i' % i)

for place in places:
    net.places.add(places[place])

transitions = {'register_request': PetriNet.Transition('a', 'register request'),
               'examine_thoroughly': PetriNet.Transition('b', 'examine thoroughly'),
               'examine_causally': PetriNet.Transition('c', 'examine causally'),
               'check_ticket': PetriNet.Transition('d', 'check ticket'),
               'decide': PetriNet.Transition('e', 'decide'),
               'reinitiate_request': PetriNet.Transition('f', 'reinitiate request'),
               'pay_compensation': PetriNet.Transition('g', 'pay compensation'),
               'reject_request': PetriNet.Transition('h', 'reject request')}

for transition in transitions:
    # print(transition)
    net.transitions.add(transitions[transition])

utils.add_arc_from_to(places['p_1'], transitions['register_request'], net)
utils.add_arc_from_to(transitions['register_request'], places['p_2'], net)
utils.add_arc_from_to(transitions['register_request'], places['p_3'], net)
utils.add_arc_from_to(places['p_2'], transitions['examine_thoroughly'], net)
utils.add_arc_from_to(places['p_2'], transitions['examine_causally'], net)
utils.add_arc_from_to(places['p_3'], transitions['check_ticket'], net)
utils.add_arc_from_to(transitions['examine_thoroughly'], places['p_4'], net)
utils.add_arc_from_to(transitions['examine_causally'], places['p_4'], net)
utils.add_arc_from_to(transitions['check_ticket'], places['p_5'], net)
utils.add_arc_from_to(places['p_4'], transitions['decide'], net)
utils.add_arc_from_to(places['p_5'], transitions['decide'], net)
utils.add_arc_from_to(transitions['decide'], places['p_6'], net)
utils.add_arc_from_to(places['p_6'], transitions['pay_compensation'], net)
utils.add_arc_from_to(places['p_6'], transitions['reject_request'], net)
utils.add_arc_from_to(places['p_6'], transitions['reinitiate_request'], net)
utils.add_arc_from_to(transitions['reinitiate_request'], places['p_2'], net)
utils.add_arc_from_to(transitions['reinitiate_request'], places['p_3'], net)
utils.add_arc_from_to(transitions['pay_compensation'], places['p_7'], net)
utils.add_arc_from_to(transitions['reject_request'], places['p_7'], net)

initial_marking = Marking()
initial_marking[places['p_1']] = 1
final_marking = Marking()
final_marking[places['p_7']] = 1

parameters = {"format": "svg", 'debug': True}
# parameters = None
gviz = pn_vis_factory.apply(net, initial_marking, final_marking, parameters=parameters)
pn_vis_factory.view(gviz)

# LOG
log = xes_import_factory.apply(os.path.join("..", "tests", "input_data", "trace_from_paper.xes"))
trace = log[0]  # abefbh
print(trace)



# Alignment computation
model_cost_function = {}
sync_cost_function = {}
for t in net.transitions:
    if t.label is not None:
        model_cost_function[t] = 1
        sync_cost_function[t] = 0
    else:
        model_cost_function[t] = 1

# giving each event in the trace the same cost
trace_costs = list(map(lambda e: 1, trace))

params = {PARAM_MODEL_COST_FUNCTION: model_cost_function, PARAM_TRACE_COST_FUNCTION: trace_costs,
          PARAM_SYNC_COST_FUNCTION: sync_cost_function}
alignment = alignments.factory.apply(trace, net, initial_marking, final_marking, parameters=params)
print("alignment length: " + str(len(alignment['alignment'])))
print(alignment)
