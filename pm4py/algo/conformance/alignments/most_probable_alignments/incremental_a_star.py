from cvxopt import matrix, solvers
from cvxopt.modeling import op, variable, constraint, sum
import numpy as np

from pm4py.algo.conformance.alignments.most_probable_alignments.probability_computation_utils import \
    calculate_model_move_probabilities_without_prior, calculate_log_move_probability, get_move_cost, \
    get_log_move_probability
from pm4py.objects import petri
from pm4py.objects.log.importer.xes.factory import import_log_from_string
from pm4py.objects.petri.petrinet import PetriNet, Marking
from pm4py.objects.petri import utils as petri_net_utils
from pm4py.visualization.petrinet import factory as petri_net_visualization_factory
from pm4py.objects.log.importer.xes import factory as xes_importer
from pm4py.util.create_artificial_event_log import create_xes_string
from pm4py.algo.conformance.alignments.utils import SKIP
from pm4py.algo.conformance.alignments.most_probable_alignments.miscellaneous_utils import apply_log_transformation, \
    is_model_move, is_log_move

model_move_probabilities_for_heuristic = None


def __compute_heuristic(sync_net, current_marking, log_move_probabilities, model_move_probabilities,
                        sync_net_final_marking):
    """
    Computes an exact heuristic using an LP based on the marking equation.

    Parameters
    ----------
    :param sync_net: synchronous product net
    :param incidence_matrix: incidence matrix
    :param current_marking: marking to start from
    :param model_move_probabilities: cost vector
    :param sync_net_final_marking: marking to reach

    Returns
    -------
    :return: h: heuristic value, x: solution vector
    """
    variables = {}
    constraints_related_to_trace_net_places = []
    constraints_related_to_process_net_places = []
    non_negativity_constraints = []
    for t in sync_net.transitions:
        # print(t.name)
        variables[t] = variable(1, str(t.name))
        # all variables must be greater or equal 0
        c = (variables[t] >= 0)
        non_negativity_constraints.append(c)
    # print("----")

    for p in sync_net.places:
        arcs_to_transitions = []  # list of all transitions that have an incoming arc from the current place
        arcs_from_transitions = []  # list of all transitions that have an arc pointing to the current place

        for out_arc in p.out_arcs:
            arcs_to_transitions.append(out_arc.target)

        for in_arc in p.in_arcs:
            arcs_from_transitions.append(in_arc.source)

        if p.name[1] == SKIP:
            print("trace net place")
            # place belongs to the trace net part
            c = (current_marking[p] + sum((-1) * variables[x] for x in arcs_to_transitions) + sum(
                variables[x] for x in arcs_from_transitions) == sync_net_final_marking[p])
            constraints_related_to_trace_net_places.append(c)
        else:
            # place belongs to the process net part
            print("process net place")
            c = (current_marking[p] + sum((-1) * variables[x] for x in arcs_to_transitions) +
                 sum(1 * variables[x] for x in arcs_from_transitions) >= 0)
            constraints_related_to_process_net_places.append(c)

        # TODO model the constraint that at least 1 token is in the process net part

        print("Outgoing arcs:")
        print(arcs_to_transitions)
        print("Incoming arcs:")
        print(arcs_from_transitions)
        print("------------------")
    # calculate the costs for each transition
    # TODO how to deal with infinite cost
    costs = {}
    for t in sync_net.transitions:
        c = get_move_cost_for_heuristic(t, log_move_probabilities, model_move_probabilities)
        costs[t] = c

    objective_function = sum(costs[x] * variables[x] for x in variables)
    print(objective_function)
    print(len(variables))
    lp = op(objective_function,
            constraints_related_to_process_net_places + constraints_related_to_trace_net_places + non_negativity_constraints)

    # TODO only for debugging issues

    print("\nConstriants")
    for c in lp.constraints():
        pass
        # print(c)
    # END
    print("\nObjective:")
    print(lp.objective)

    equalities = lp.equalities()
    inequalities = lp.inequalities()
    print(equalities)
    print(inequalities)

    lp.solve()
    print(lp.status)
    print("-----------------")
    print("Variables")
    print(lp.variables())
    for v in variables:
        print(variables[v].value)

    print("Objective value: " + str(lp.objective.value()))
    return


def get_move_cost_for_heuristic(t, log_move_probabilities, model_move_probabilities):
    if is_model_move(t, SKIP):
        return apply_log_transformation(get_model_move_probability_ignoring_marking(t, model_move_probabilities))
    elif is_log_move(t, SKIP):
        return apply_log_transformation(get_log_move_probability(t, log_move_probabilities))
    else:
        cost = min(
            apply_log_transformation(get_log_move_probability(t, log_move_probabilities)),
            apply_log_transformation(
                get_model_move_probability_ignoring_marking(t, model_move_probabilities)))
        return cost


def get_model_move_probability_ignoring_marking(requested_transition, model_move_probabilities):
    """
    the dict model_move_probabilities holds conditional probabilities of the form P(T=t|M=m), i.e., given a marking how
    likely is it that transition t was executed. Since we cannot handle these conditional probabilities/ costs in the
    heuristic calculation, this function computes the sample mean of the transitions ignoring the marking.

    :param requested_transition: Transition object
    :param model_move_probabilities: dict with markings as keys and transitions and their frequencies/ probabilities as
                                     value
    :return: probability (ignoring the conditional dependence) that the requested_transition was executed
    """
    global model_move_probabilities_for_heuristic

    if not model_move_probabilities_for_heuristic:
        model_move_probabilities_for_heuristic = {}
        sum_transitions_frequencies = 0
        for marking in model_move_probabilities:
            for t in marking["outgoing_transitions"]:
                # t is a dict like:
                # <class 'dict'>: {'unique_name': 't_A', 'frequency': 50.0, 'model_move_probability': 1.0}
                if t["unique_name"] in model_move_probabilities_for_heuristic:
                    model_move_probabilities_for_heuristic[t["unique_name"]]["frequency"] += \
                        t["frequency"]
                else:
                    model_move_probabilities_for_heuristic[t["unique_name"]] = {
                        "frequency": t["frequency"]}
                sum_transitions_frequencies += t["frequency"]
        for t in model_move_probabilities_for_heuristic:
            model_move_probabilities_for_heuristic[t]["probability"] = model_move_probabilities_for_heuristic[t][
                                                                           "frequency"] / sum_transitions_frequencies

    # TODO remove
    debugging_model_move_probabilities_for_heuristic = model_move_probabilities_for_heuristic
    # get the corresponding transition probability
    keys = model_move_probabilities_for_heuristic.keys()
    if requested_transition.name[1] in model_move_probabilities_for_heuristic:
        return model_move_probabilities_for_heuristic[requested_transition.name[1]]["probability"]
    else:
        # requested_transition was never executed
        return 0


if __name__ == '__main__':
    print(sum(1 * x for x in [2, 5, 7, 8]))
    print(sum([]))
    # TEST ####### TEST ####### TEST ####### TEST ####### TEST ####### TEST ####### TEST ####### TEST ####### TEST #####

    # create petri net
    test_petri_net = PetriNet("test")
    places = {}
    for i in range(1, 4):
        places['p_%i' % i] = PetriNet.Place('p_%i' % i)
        test_petri_net.places.add(places['p_%i' % i])

    transitions = {
        # PetriNet.Transition(<unique name in petri net>, <label>)
        't_A': PetriNet.Transition('t_A', 'A'),
        't_B': PetriNet.Transition('t_B', 'B'),
        't_C': PetriNet.Transition('t_C', 'C'),
        # 't_D': PetriNet.Transition('t_D', 'None'),
        # 't_E': PetriNet.Transition('t_E', 'None')

    }

    for transition in transitions:
        test_petri_net.transitions.add(transitions[transition])

    petri_net_utils.add_arc_from_to(places['p_1'], transitions['t_A'], test_petri_net)
    petri_net_utils.add_arc_from_to(transitions['t_A'], places['p_2'], test_petri_net)
    petri_net_utils.add_arc_from_to(places['p_2'], transitions['t_B'], test_petri_net)
    petri_net_utils.add_arc_from_to(places['p_2'], transitions['t_C'], test_petri_net)
    petri_net_utils.add_arc_from_to(transitions['t_B'], places['p_3'], test_petri_net)
    petri_net_utils.add_arc_from_to(transitions['t_C'], places['p_3'], test_petri_net)

    initial_marking = Marking()
    initial_marking[places['p_1']] = 1
    final_marking = Marking()
    final_marking[places['p_3']] = 1

    traces = [
        {"frequency": 1, "events": ["A"]}
    ]
    event_log = xes_importer.import_log_from_string(create_xes_string(traces))

    trace_net, trace_im, trace_fm = petri.utils.construct_trace_net(event_log[0])
    sync_prod_net, sync_initial_marking, sync_final_marking = petri.synchronous_product.construct(trace_net, trace_im,
                                                                                                  trace_fm,
                                                                                                  test_petri_net,
                                                                                                  initial_marking,
                                                                                                  final_marking,
                                                                                                  SKIP)

    # gviz = petri_net_visualization_factory.apply(sync_prod_net, sync_initial_marking, sync_final_marking,
    #                                             parameters={"format": "svg", 'debug': True})
    # petri_net_visualization_factory.view(gviz)
    # incidence_matrix = petri.incidence_matrix.construct(sync_prod_net)

    traces = [
        {"frequency": 40, "events": ["A", "B"]},
        {"frequency": 10, "events": ["A", "C"]}
    ]
    xes_str = create_xes_string(traces)
    event_log = import_log_from_string(xes_str)

    log_move_prior = None
    log_move_prob = calculate_log_move_probability(event_log, log_move_prior)

    model_move_probabilities_without_prior = calculate_model_move_probabilities_without_prior(event_log, test_petri_net,
                                                                                              initial_marking,
                                                                                              final_marking)

    __compute_heuristic(sync_prod_net, sync_initial_marking, log_move_prob,
                        model_move_probabilities_without_prior, sync_final_marking)
