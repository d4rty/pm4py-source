import os

from pm4py import util
from pm4py.algo.conformance import alignments as ali
from pm4py.algo.conformance.alignments.versions.state_equation_a_star import PARAM_MODEL_COST_FUNCTION
from pm4py.algo.conformance.alignments.versions.state_equation_a_star import PARAM_SYNC_COST_FUNCTION
from pm4py.algo.conformance.alignments.versions.state_equation_a_star import PARAM_TRACE_COST_FUNCTION
from pm4py.objects import log as log_lib
from pm4py.objects import petri as petri
from pm4py.objects.log.importer.xes import factory as xes_importer
from pm4py.visualization.petrinet import factory
from pm4py.objects.petri import semantics


def align(trace, net, im, fm, model_cost_function, sync_cost_function):
    trace_costs = list(map(lambda e: 1000, trace))
    params = dict()
    params[util.constants.PARAMETER_CONSTANT_ACTIVITY_KEY] = log_lib.util.xes.DEFAULT_NAME_KEY
    params[PARAM_MODEL_COST_FUNCTION] = model_cost_function
    params[PARAM_TRACE_COST_FUNCTION] = trace_costs
    params[PARAM_SYNC_COST_FUNCTION] = sync_cost_function
    return ali.factory.apply_trace(trace, net, im, fm, parameters=params,
                                   version=ali.factory.VERSION_STATE_EQUATION_A_STAR)


def execute_script():
    log_path = os.path.join("..", "tests", "input_data", "running-example.xes")
    pnml_path = os.path.join("..", "tests", "input_data", "running-example.pnml")

    # log_path = 'C:/Users/bas/Documents/tue/svn/private/logs/a32_logs/a32f0n05.xes'
    # pnml_path = 'C:/Users/bas/Documents/tue/svn/private/logs/a32_logs/a32.pnml'

    log = xes_importer.import_log(log_path)

    # for event in log:
    #    print(type(event))
    #    print(event)

    net, marking, fmarking = petri.importer.pnml.import_net(
        pnml_path)

    # shows the enabled transitions for a particular marking
    transitions = semantics.enabled_transitions(net, marking)

    # shows the unique label of the places
    # print(net.places)
    # for place in net.places:
    #     print(place.name)
    # print(net.transitions)
    # for transition in net.transitions:
    #     print("\nTRANS: ", transition.name, transition.label)

    model_cost_function = dict()
    sync_cost_function = dict()
    for t in net.transitions:
        if t.label is not None:
            model_cost_function[t] = 1000
            sync_cost_function[t] = 0
        else:
            model_cost_function[t] = 1

    # visualization of petri net
    # factory.view(factory.apply(net))
    # factory.view(factory.apply(net, marking, fmarking))

    # print(list(map(lambda trace: align(trace, net, marking, fmarking, model_cost_function, sync_cost_function), log)))

    alignment = ali.factory.apply(log[0], net, marking, fmarking)

    for move in alignment['alignment']:
        print("Label:   " + str(move['label']))
        print("Name:    " + str(move['name']))
        print("Marking: " + str(move['marking']))
        print()


if __name__ == '__main__':
    execute_script()
