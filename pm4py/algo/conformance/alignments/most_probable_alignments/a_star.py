"""
This module contains code that allows us to compute alignments on the basis of a regular A* search on the state-space
of the synchronous product net of a trace and a Petri net.
The main algorithm follows [1]_.
When running the log-based variant, the code is running in parallel on a trace based level.
Furthermore, by default, the code applies heuristic estimation, and prefers those states that have the smallest h-value
in case the f-value of two states is equal.

References
----------
.. [1] Sebastiaan J. van Zelst et al., "Tuning Alignment Computation: An Experimental Evaluation",
      ATAED@Petri Nets/ACSD 2017: 6-20. `http://ceur-ws.org/Vol-1847/paper01.pdf`_.

"""
import math
import heapq
from typing import Any

from dataclasses import dataclass

from pm4py import util as pm4pyutil
from pm4py.algo.conformance import alignments
from pm4py.objects import petri
from pm4py.objects.log.util.xes import DEFAULT_NAME_KEY
from pm4py.util.constants import PARAMETER_CONSTANT_ACTIVITY_KEY
from pm4py.visualization.petrinet import factory as pn_vis_factory

PARAM_TRACE_COST_FUNCTION = 'trace_cost_function'
PARAM_MODEL_COST_FUNCTION = 'model_cost_function'
PARAM_SYNC_COST_FUNCTION = 'sync_cost_function'

PARAMETERS = [PARAM_TRACE_COST_FUNCTION, PARAM_MODEL_COST_FUNCTION, PARAM_SYNC_COST_FUNCTION,
              pm4pyutil.constants.PARAMETER_CONSTANT_ACTIVITY_KEY]


def apply(trace, petri_net, initial_marking, final_marking, log_move_probabilities, model_move_probabilities,
          parameters=None):
    activity_key = DEFAULT_NAME_KEY if parameters is None or PARAMETER_CONSTANT_ACTIVITY_KEY not in parameters else \
        parameters[pm4pyutil.constants.PARAMETER_CONSTANT_ACTIVITY_KEY]

    trace_net, trace_initial_marking, trace_final_marking = petri.utils.construct_trace_net(trace,
                                                                                            activity_key=activity_key)

    sync_prod, sync_initial_marking, sync_final_marking = petri.synchronous_product.construct(trace_net,
                                                                                              trace_initial_marking,
                                                                                              trace_final_marking,
                                                                                              petri_net,
                                                                                              initial_marking,
                                                                                              final_marking,
                                                                                              alignments.utils.SKIP)

    # view synchronous product net
    gviz = pn_vis_factory.apply(sync_prod, sync_initial_marking, sync_final_marking,
                                parameters={"debug": True, "format": "svg"})
    # pn_vis_factory.view(gviz)

    return apply_sync_prod(sync_prod, petri_net, sync_initial_marking, sync_final_marking, log_move_probabilities,
                           model_move_probabilities, alignments.utils.SKIP)


def apply_sync_prod(sync_prod, process_net, initial_marking, final_marking, log_move_probabilities,
                    model_move_probabilities, skip):
    """
    Performs the basic alignment search on top of the synchronous product net, given a cost function and skip-symbol

    Parameters
    ----------
    sync_prod: :class:`pm4py.objects.petri.net.PetriNet` synchronous product net
    initial_marking: :class:`pm4py.objects.petri.net.Marking` initial marking in the synchronous product net
    final_marking: :class:`pm4py.objects.petri.net.Marking` final marking in the synchronous product net
    skip: :class:`Any` symbol to use for skips in the alignment

    Returns
    -------
    dictionary : :class:`dict` with keys **alignment**, **cost**, **visited_states**, **queued_states**
    and **traversed_arcs**
    """
    return __search(sync_prod, process_net, initial_marking, final_marking, log_move_probabilities,
                    model_move_probabilities, skip)


def __search(sync_net, process_net, initial_marking, final_marking, log_move_probabilities, model_move_probabilities,
             skip):
    incidence_matrix = petri.incidence_matrix.construct(sync_net)

    closed = set()
    h = __compute_exact_heuristic(sync_net, incidence_matrix, initial_marking)
    ini_state = SearchTuple(0 + h, 0, h, initial_marking, None, None, 1)
    open_set = [ini_state]  # visited markings
    visited = 0
    queued = 0
    traversed = 0
    while not len(open_set) == 0:
        current_state = heapq.heappop(open_set)
        visited += 1
        current_marking = current_state.m
        closed.add(current_marking)

        if current_marking == final_marking:
            return __reconstruct_alignment(current_state, visited, queued, traversed)

        for t in petri.semantics.enabled_transitions(sync_net, current_marking):
            if current_state.t is not None and __is_log_move(current_state.t, skip) and __is_model_move(t, skip):
                # TODO why?? maybe the node sequence (LOG_MOVE->MODEL_MOVE) can't be optimal, therefore skip it
                continue
            traversed += 1
            new_marking = petri.semantics.execute(t, sync_net, current_marking)
            if new_marking in closed:
                continue

            cost = __get_move_cost(t, current_marking, log_move_probabilities, model_move_probabilities, process_net)
            g = current_state.g + cost

            alt = next((enum[1] for enum in enumerate(open_set) if enum[1].m == new_marking), None)
            if alt is not None:
                if g >= alt.g:
                    continue
                open_set.remove(alt)
                heapq.heapify(open_set)
            queued += 1
            probability = current_state.probability * __get_move_probability(t, current_marking, log_move_probabilities,
                                                                             model_move_probabilities)
            h, x = __derive_heuristic(incidence_matrix, None, None, t, current_state.h)
            tp = SearchTuple(g + h, g, h, new_marking, current_state, t, probability)

            heapq.heappush(open_set, tp)
            heapq.heapify(open_set)


def __reconstruct_alignment(state, visited, queued, traversed):
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
    return {'alignment': alignment, 'cost': state.g, 'probability': math.exp(state.g), 'visited_states': visited,
            'queued_states': queued,
            'traversed_arcs': traversed, "probability": state.probability}


def __derive_heuristic(incidence_matrix, cost_vec, x, t, h):
    return max(0, 0), 1


def __is_model_move(t, skip):
    return t.label[0] == skip and t.label[1] != skip


def __is_log_move(t, skip):
    return t.label[0] != skip and t.label[1] == skip


def __get_move_cost(transition, marking, log_move_probabilities, model_move_probabilities, process_net):
    if __is_model_move(transition, alignments.utils.SKIP):
        return __apply_log_transformation(
            __get_model_move_probability(transition, marking, model_move_probabilities, process_net))
    elif __is_log_move(transition, alignments.utils.SKIP):
        return __apply_log_transformation(__get_log_move_probability(transition, log_move_probabilities))
    else:
        # synchronous move
        cost = min(
            __apply_log_transformation(__get_log_move_probability(transition, log_move_probabilities)),
            __apply_log_transformation(
                __get_model_move_probability(transition, marking, model_move_probabilities, process_net)))
        return cost


def __get_move_probability(transition, marking, log_move_probabilities, model_move_probabilities, process_net):
    if __is_model_move(transition, alignments.utils.SKIP):
        return __get_model_move_probability(transition, marking, model_move_probabilities, process_net)
    elif __is_log_move(transition, alignments.utils.SKIP):
        return __get_log_move_probability(transition, log_move_probabilities, process_net)
    else:
        # synchronous move
        cost = max(__get_log_move_probability(transition, log_move_probabilities),
                   __get_model_move_probability(transition, marking, model_move_probabilities, process_net))
        return cost


def __get_model_move_probability(transition, marking_sync_product_net, model_move_probabilities, process_net):
    def is_model_place(p):
        return p.name[0] == alignments.utils.SKIP and p.name[1] != alignments.utils.SKIP

    marking_of_process_net = {}
    for place in marking_sync_product_net:
        if is_model_place(place):
            number_tokens = marking_sync_product_net[place]
            process_net_place_name = place.name[1]
            marking_of_process_net[process_net_place_name] = number_tokens

    model_move_probabilities_for_marking = [d for d in model_move_probabilities if
                                            d['marking'] == marking_of_process_net]

    if len(model_move_probabilities_for_marking) == 1:
        # current marking exists in model_move_probabilities
        for tr in model_move_probabilities_for_marking[0]['outgoing_transitions']:
            if tr['unique_name'] == transition.name[1]:
                return tr['model_move_probability']
        # no probability for given transition
        # TODO apply prior?!
        return 0

    elif len(model_move_probabilities_for_marking) == 0:
        # current marking does not exist in model_move_probabilities
        # TODO apply prior?!
        return 0
    else:
        # something went wrong, markings should be unique
        raise Exception('Multiple markings founded')


def calculate_model_move_prior():
    pass


def __get_log_move_probability(transition, log_move_probabilities):
    return log_move_probabilities[transition.label[0]]


def __compute_exact_heuristic(sync_net, incidence_matrix, marking):
    # TODO
    # smallest floating-point value after 0
    return 0


def __apply_log_transformation(probability):
    # TODO optimize --> apply log transformation not on the fly but apply it if possible before
    # --> try to avoid duplicate calls
    if not (0 <= probability <= 1):
        raise ValueError('probability needs to be in the range [0,1]')
    if probability > 0:
        return -math.log(probability)
    elif probability == 0:
        # -log(0) ~ math.inf
        return math.inf
    else:
        raise ValueError('Log Transformation not applicable to numbers below 0')


@dataclass
class SearchTuple:
    f: float  # f = g + h value corresponding to a* algorithm
    g: float  # g (costs to this node) value corresponding to a* algorithm
    h: float  # h (estimated costs from this node to target node) value corresponding to a* algorithm
    m: petri.petrinet.Marking
    p: Any  # previous search tuple/ parent
    t: petri.petrinet.PetriNet.Transition
    probability: float  # probability of reaching this state

    # x < y --> calls x.__lt__(y)
    def __lt__(self, other):
        # returns true if x.f is > y.f
        # (reason: always have the SearchTuple with highest f value on heapq at first position)
        if self.f < other.f:
            return True
        elif other.f < self.f:
            return False
        else:
            if self.h < other.h:
                return True

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
