import math
import heapq
from typing import Any
from dataclasses import dataclass
from pm4py import util as pm4pyutil
from pm4py.algo.conformance import alignments
from pm4py.algo.conformance.alignments.most_probable_alignments.probability_computation_utils import get_move_cost, \
    get_move_probability
from pm4py.objects import petri
from pm4py.objects.log.util.xes import DEFAULT_NAME_KEY
from pm4py.util.constants import PARAMETER_CONSTANT_ACTIVITY_KEY
from pm4py.visualization.petrinet import factory as pn_vis_factory
from pm4py.algo.conformance.alignments.most_probable_alignments.miscellaneous_utils import is_log_move, is_model_move


def apply(trace, petri_net, initial_marking, final_marking, log_move_probabilities, model_move_probabilities,
          parameters=None):
    print(trace)
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

    gviz = pn_vis_factory.apply(sync_prod,
                                sync_initial_marking,
                                sync_final_marking,
                                parameters={"debug": True, "format": "svg"})
    # pn_vis_factory.view(gviz)

    return apply_sync_prod(sync_prod, petri_net, sync_initial_marking, sync_final_marking, log_move_probabilities,
                           model_move_probabilities, alignments.utils.SKIP)


def apply_sync_prod(sync_prod, process_net, initial_marking, final_marking, log_move_probabilities,
                    model_move_probabilities, skip):
    """

    :param sync_prod:
    :param process_net:
    :param initial_marking:
    :param final_marking:
    :param log_move_probabilities:
    :param model_move_probabilities:
    :param skip:
    :return:
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
            if current_state.t is not None and is_log_move(current_state.t, skip) and is_model_move(t, skip):
                # TODO why?? maybe the node sequence (LOG_MOVE->MODEL_MOVE) can't be optimal, therefore skip it?!
                continue
            traversed += 1
            new_marking = petri.semantics.execute(t, sync_net, current_marking)
            if new_marking in closed:
                continue

            cost = get_move_cost(t, current_marking, log_move_probabilities, model_move_probabilities, process_net)
            g = current_state.g + cost

            alt = next((enum[1] for enum in enumerate(open_set) if enum[1].m == new_marking), None)
            if alt is not None:
                if g >= alt.g:
                    continue
                open_set.remove(alt)
                heapq.heapify(open_set)
            queued += 1
            probability = current_state.probability * get_move_probability(t, current_marking, log_move_probabilities,
                                                                           model_move_probabilities, process_net)
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


def __compute_exact_heuristic(sync_net, incidence_matrix, marking):
    return 0


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
