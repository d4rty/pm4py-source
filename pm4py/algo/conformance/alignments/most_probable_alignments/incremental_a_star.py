import heapq

import math

from ortools.linear_solver import pywraplp

import pm4py
from pm4py.algo.conformance import alignments
from pm4py.algo.conformance.alignments.most_probable_alignments.probability_computation_utils import \
    calculate_model_move_probabilities_without_prior, calculate_log_move_probability, get_move_cost, \
    get_log_move_probability
from pm4py.evaluation.replay_fitness.versions.alignment_based import DEFAULT_NAME_KEY
from pm4py.objects import petri
from pm4py.objects.log.importer.xes.factory import import_log_from_string
from pm4py.objects.log.log import Trace
from pm4py.objects.petri.petrinet import PetriNet, Marking
from pm4py.objects.petri import utils as petri_net_utils
from pm4py.objects.petri.utils import construct_trace_net_cost_aware
from pm4py.util.constants import PARAMETER_CONSTANT_ACTIVITY_KEY
from pm4py.visualization.petrinet import factory as petri_net_visualization_factory
from pm4py.objects.log.importer.xes import factory as xes_importer
from pm4py.util.create_artificial_event_log import create_xes_string
from pm4py.algo.conformance.alignments.utils import SKIP, print_alignment
from pm4py.algo.conformance.alignments.most_probable_alignments.miscellaneous_utils import apply_log_transformation, \
    is_model_move, is_log_move, place_from_synchronous_product_net_belongs_to_process_net_part, \
    place_from_synchronous_product_net_belongs_to_trace_net_part
from pm4py.visualization.petrinet import factory as pn_vis_factory

model_move_probabilities_for_heuristic = None


def apply(trace, petri_net, initial_marking, final_marking, parameters=None):
    activity_key = DEFAULT_NAME_KEY if parameters is None or PARAMETER_CONSTANT_ACTIVITY_KEY not in parameters else \
        parameters[
            pm4py.util.constants.PARAMETER_CONSTANT_ACTIVITY_KEY]
    incremental_trace = Trace()

    # create empty closed and open set
    open_set = []
    closed_set = set()
    first_event = True
    for event in trace:
        incremental_trace.append(event)
        print(incremental_trace)
        if first_event:
            print("first event")
            # activity_key: :class:`str` key of the attribute of the events that defines the activity name
            trace_net, trace_im, trace_fm = petri.utils.construct_trace_net(incremental_trace,
                                                                            activity_key=activity_key)
            sync_prod, sync_im, sync_fm = petri.synchronous_product.construct(trace_net,
                                                                              trace_im,
                                                                              trace_fm,
                                                                              petri_net,
                                                                              initial_marking,
                                                                              final_marking,
                                                                              alignments.utils.SKIP)
            first_event = False
        else:
            print("not first event")
            sync_prod, sync_fm = petri.synchronous_product.extend_trace_net_of_synchronous_product_net(sync_prod, event,
                                                                                                       sync_fm, SKIP,
                                                                                                       activity_key)

        gviz = pn_vis_factory.apply(sync_prod, sync_im, sync_fm,
                                    parameters={"debug": True, "format": "svg"})
        pn_vis_factory.view(gviz)

        cost_function = alignments.utils.construct_standard_cost_function(sync_prod, alignments.utils.SKIP)

        prefix_alignment, open_set, closed_set = __search(sync_prod, sync_im, sync_fm,
                                                          cost_function,
                                                          alignments.utils.SKIP, open_set, closed_set)
        print(prefix_alignment)
        print_alignment(prefix_alignment)
        print(open_set)
        print(closed_set)


def __search(sync_net, ini, fin, cost_function, skip, open_set_heap, closed_set):
    print("hello")
    # TODO
    h, x = 0, True
    ini_state = SearchTuple(0 + h, 0, h, ini, None, None, x, True)
    if len(open_set_heap) == 0:
        open_set_heap = [ini_state]
    else:
        # TODO recalculate heurisitc!!
        heapq.heapify(open_set_heap)  # visited markings
    visited = 0
    queued = 0
    traversed = 0
    while not len(open_set_heap) == 0:
        curr = heapq.heappop(open_set_heap)
        if not curr.trust:
            h, x = 0, True
            tp = SearchTuple(curr.g + h, curr.g, h, curr.m, curr.p, curr.t, x, True)
            heapq.heappush(open_set_heap, tp)
            heapq.heapify(open_set_heap)  # transform a populated list into a heap
            continue

        visited += 1
        current_marking = curr.m
        closed_set.add(current_marking)

        # check if we reached the final marking

        for place in current_marking:
            if place_from_synchronous_product_net_belongs_to_trace_net_part(place):
                for place2 in fin:
                    if place_from_synchronous_product_net_belongs_to_trace_net_part(place2):
                        if place.name == place2.name:
                            return __reconstruct_alignment(curr, visited, queued, traversed), open_set_heap, closed_set

        print("enabled transitions:")
        for t in petri.semantics.enabled_transitions(sync_net, current_marking):
            print(t)
            if curr.t is not None and is_log_move(curr.t, skip) and is_model_move(t, skip):
                continue

            traversed += 1
            new_marking = petri.semantics.execute(t, sync_net, current_marking)
            if new_marking in closed_set:
                continue
            g = curr.g + cost_function[t]

            # enum is a tuple (int, SearchTuple), alt is a SearchTuple
            alt = next((enum[1] for enum in enumerate(open_set_heap) if enum[1].m == new_marking), None)
            if alt is not None:
                if g >= alt.g:
                    continue
                open_set_heap.remove(alt)
                heapq.heapify(open_set_heap)
            queued += 1
            h, x = 0, True
            tp = SearchTuple(g + h, g, h, new_marking, curr, t, x, True)
            heapq.heappush(open_set_heap, tp)
            heapq.heapify(open_set_heap)


def __reconstruct_alignment(state, visited, queued, traversed):
    # state is a SearchTuple
    parent = state.p
    alignment = [{"marking_before_transition": state.p.m,
                  "label": state.t.label,
                  "name": state.t.name,
                  "marking_after_transition": state.m}]
    while parent.p is not None:
        alignment = [{"marking_before_transition": parent.p.m,
                      "label": parent.t.label,
                      "name": parent.t.name,
                      "marking_after_transition": parent.m}] + alignment
        parent = parent.p
    return {'alignment': alignment, 'cost': state.g, 'visited_states': visited, 'queued_states': queued,
            'traversed_arcs': traversed}


def __compute_heuristic_most_probable_alignments(sync_net, current_marking, log_move_probabilities,
                                                 model_move_probabilities,
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
    costs = {}
    for t in sync_net.transitions:
        c = get_move_cost_for_heuristic(t, log_move_probabilities, model_move_probabilities)
        costs[t] = c
    solver = pywraplp.Solver('SolveSimpleSystem', pywraplp.Solver.GLOP_LINEAR_PROGRAMMING)
    variables = {}
    constraints = []
    for t in sync_net.transitions:
        if costs[t] < math.inf:
            # only create variables that have finite cost/ probability > 0
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
    # store coefficients for each variable here because when calling constraint.SetCoefficient multiple times for the
    # same variable it overwrites always the last value for the given variable, i.e. it is NOT possible to model the
    # following constraint: x >= x1 + x2 -x1 with:
    # c.SetCoefficient(x1 , 1)
    # c.SetCoefficient(x2 , 1)
    # c.SetCoefficient(x1 , -1) --> overwrites the previous coefficient of x1
    constraint_one_token_in_process_net_part_coefficients = {}
    for v in variables:
        constraint_one_token_in_process_net_part_coefficients[v] = 0

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
            c = solver.Constraint(lb_and_ub, lb_and_ub)
        else:
            # place belongs to the process net part
            # enforce that the constraint is greater or equal 0, i.e.,
            # constraint + constant >= 0  -->  constraint >= 0 - constant
            c = solver.Constraint(0 - current_marking[p], solver.infinity())

            for t in arcs_to_transitions:
                if t in variables:
                    constraint_one_token_in_process_net_part_coefficients[t] -= 1

            for t in arcs_from_transitions:
                if t in variables:
                    constraint_one_token_in_process_net_part_coefficients[t] += 1

        for t in arcs_to_transitions:
            if t in variables:
                c.SetCoefficient(variables[t], -1)
        for t in arcs_from_transitions:
            if t in variables:
                c.SetCoefficient(variables[t], 1)
        constraints.append(c)
    # build constraint that enforces at least one token in the process net part
    for v in variables:
        constraint_one_token_in_process_net_part.SetCoefficient(variables[v],
                                                                constraint_one_token_in_process_net_part_coefficients[v]
                                                                )
    objective = solver.Objective()
    for v in variables:
        objective.SetCoefficient(variables[v], costs[v])
    objective.SetMinimization()
    solver.Solve()
    # debugging
    # print('Number of variables =', solver.NumVariables())
    # print('Number of constraints =', solver.NumConstraints())
    # print('Solution:')
    # for v in variables:
    #     print(str(v.name) + ":" + str(variables[v].solution_value()))
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


class SearchTuple:
    def __init__(self, f, g, h, m, p, t, x, trust):
        self.f = f
        self.g = g
        self.h = h
        self.m = m
        self.p = p
        self.t = t
        self.x = x
        self.trust = trust

    def __lt__(self, other):
        if self.f < other.f:
            return True
        elif other.f < self.f:
            return False
        else:
            return self.h < other.h

    def __get_firing_sequence(self):
        ret = []
        if self.p is not None:
            ret = ret + self.p.__get_firing_sequence()
        if self.t is not None:
            ret.append(self.t)
        return ret

    def __repr__(self):
        string_build = ["\nm=" + str(self.m), " f=" + str(self.f), ' g=' + str(self.g), " h=" + str(self.h),
                        " path=" + str(self.__get_firing_sequence()) + "\n\n"]
        return " ".join(string_build)


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

    gviz = petri_net_visualization_factory.apply(test_petri_net, sync_initial_marking, sync_final_marking,
                                                 parameters={"format": "svg", 'debug': True})
    petri_net_visualization_factory.view(gviz)
    incidence_matrix = petri.incidence_matrix.construct(sync_prod_net)

    traces = [
        {"frequency": 40, "events": ["A","B","C"]},
        {"frequency": 10, "events": ["A", "C"]}
    ]

    xes_str = create_xes_string(traces)
    event_log = import_log_from_string(xes_str)

    log_move_prior = None
    log_move_prob = calculate_log_move_probability(event_log, log_move_prior)

    model_move_probabilities_without_prior = calculate_model_move_probabilities_without_prior(event_log, test_petri_net,
                                                                                              initial_marking,
                                                                                              final_marking)

    res = __compute_heuristic_most_probable_alignments(sync_prod_net, sync_initial_marking, log_move_prob,
                                                       model_move_probabilities_without_prior, sync_final_marking)
    print(res)

    apply(event_log[0], test_petri_net, initial_marking, final_marking)
