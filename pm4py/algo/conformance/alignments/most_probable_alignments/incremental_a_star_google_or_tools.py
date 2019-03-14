from cvxopt import matrix, solvers
from ortools.linear_solver import pywraplp
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

    solver = pywraplp.Solver('SolveSimpleSystem', pywraplp.Solver.GLOP_LINEAR_PROGRAMMING)

    variables = {}
    constraints = []

    for t in sync_net.transitions:
        # print(t.name)
        variables[t] = solver.NumVar(0, solver.infinity(), str(t.name))

    # calculate current number of tokens in the process net part of the synchronous product net
    number_tokens_in_process_net_part = 0
    for p in current_marking:
        if p.name[0] == SKIP:
            number_tokens_in_process_net_part += current_marking[p]

    # constraint that enforces that at least one token is in the process net part of the synchronous product net
    # example: 1 <= var1 * coefficient1 + var2 * coefficient2 + ... + constant
    # rewrite to -->  1 - constant <= var1 * coefficient1 + var2 * coefficient2 + ...
    lb = 1 - number_tokens_in_process_net_part
    constraint_one_token_in_process_net_part = solver.Constraint(lb, solver.infinity())

    # define constraints
    for p in sync_net.places:
        arcs_to_transitions = []  # list of all transitions that have an incoming arc from the current place
        arcs_from_transitions = []  # list of all transitions that have an arc pointing to the current place

        for out_arc in p.out_arcs:
            arcs_to_transitions.append(out_arc.target)

        for in_arc in p.in_arcs:
            arcs_from_transitions.append(in_arc.source)

        if p.name[1] == SKIP:
            # place belongs to the trace net part
            lb_and_ub = sync_net_final_marking[p] - current_marking[p]
            # enforce that the constraint is equal to the value of lb_and_ub, i.e., constraint = lb_and_ub
            c = solver.Constraint(lb_and_ub, lb_and_ub)
        else:
            # place belongs to the process net part
            # enforce that the constraint is greater or equal 0, i.e.,
            # constraint + constant >= 0 --> constraint >= 0 - constant
            c = solver.Constraint(0 - current_marking[p], solver.infinity())

            for x in arcs_to_transitions:
                constraint_one_token_in_process_net_part.SetCoefficient(variables[x], -1)
            for x in arcs_from_transitions:
                constraint_one_token_in_process_net_part.SetCoefficient(variables[x], 1)

        for x in arcs_to_transitions:
            c.SetCoefficient(variables[x], -1)
        for x in arcs_from_transitions:
            c.SetCoefficient(variables[x], 1)
        constraints.append(c)

    # calculate the costs for each transition
    # TODO how to deal with infinite cost?!
    #  -> only add variable if costs below infinity otherwise variable cannot be part of the solution

    costs = {}
    for t in sync_net.transitions:
        c = get_move_cost_for_heuristic(t, log_move_probabilities, model_move_probabilities)
        costs[t] = c

    objective = solver.Objective()

    #    objective_function = sum(costs[x] * variables[x] for x in variables)
    for x in variables:
        objective.SetCoefficient(variables[x], costs[x])
    objective.SetMinimization()
    solver.Solve()

    print('Number of variables =', solver.NumVariables())
    print('Number of constraints =', solver.NumConstraints())
    print('Solution:')
    for v in variables:
        print(str(v.name) + ":" + str(variables[v].solution_value()))

    lp_solution = 0
    for v in variables:
        lp_solution += variables[v].solution_value() * costs[v]
    return lp_solution


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

    if requested_transition.name[1] in model_move_probabilities_for_heuristic:
        return model_move_probabilities_for_heuristic[requested_transition.name[1]]["probability"]
    else:
        # requested_transition was never executed
        return 0


if __name__ == '__main__':
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

    print(__compute_heuristic(sync_prod_net, sync_initial_marking, log_move_prob,
                        model_move_probabilities_without_prior, sync_final_marking))
