import math
import heapq
from typing import Any

from dataclasses import dataclass

from pm4py import util as pm4pyutil
from pm4py.algo.conformance import alignments
from pm4py.objects import petri
from pm4py.objects.log.util.xes import DEFAULT_NAME_KEY
from pm4py.objects.petri.petrinet import Marking
from pm4py.util.constants import PARAMETER_CONSTANT_ACTIVITY_KEY
from pm4py.visualization.petrinet import factory as pn_vis_factory


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

    # view synchronous product net
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
            if current_state.t is not None and __is_log_move(current_state.t, skip) and __is_model_move(t, skip):
                # TODO why?? maybe the node sequence (LOG_MOVE->MODEL_MOVE) can't be optimal, therefore skip it?!
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


def __is_model_move(t, skip):
    return t.label[0] == skip and t.label[1] != skip


def __is_log_move(t, skip):
    return t.label[0] != skip and t.label[1] == skip


def __get_move_cost(transition, marking, log_move_probabilities, model_move_probabilities, process_net):
    """

    :param transition:
    :param marking:
    :param log_move_probabilities:
    :param model_move_probabilities:
    :param process_net:
    :return:
    """
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
    """

    :param transition:
    :param marking:
    :param log_move_probabilities:
    :param model_move_probabilities:
    :param process_net:
    :return:
    """
    if __is_model_move(transition, alignments.utils.SKIP):
        return __get_model_move_probability(transition, marking, model_move_probabilities, process_net)
    elif __is_log_move(transition, alignments.utils.SKIP):
        return __get_log_move_probability(transition, log_move_probabilities)
    else:
        # synchronous move
        cost = max(__get_log_move_probability(transition, log_move_probabilities),
                   __get_model_move_probability(transition, marking, model_move_probabilities, process_net))
        return cost


def __sync_process_net_place_is_model_place(p):
    """
    :param p: Place object; represents a place of a synchronous product net
    :return: Boolean - true if the synchronous product net state represents a state of the process net
    """
    return p.name[0] == alignments.utils.SKIP and p.name[1] != alignments.utils.SKIP


def __get_model_move_probability(transition, marking_sync_product_net, model_move_probabilities_without_prior,
                                 process_net):
    """

    :param transition: Transition object
    :param marking_sync_product_net:
    :param model_move_probabilities_without_prior:
    :param process_net:
    :return:
    """
    constant_prior = 1
    wanted_transition = transition.name[1]

    marking_process_net = Marking()
    marking_process_net_dict = {}  # dict of marking - key: place name, value: number tokens

    for place in marking_sync_product_net:
        if __sync_process_net_place_is_model_place(place):
            number_tokens = marking_sync_product_net[place]
            process_net_place_name = place.name[1]
            # create marking object of process net
            for p in process_net.places:
                if p.name == process_net_place_name:
                    marking_process_net[p] = number_tokens
            marking_process_net_dict[process_net_place_name] = number_tokens

    model_move_probabilities_for_marking = [d for d in model_move_probabilities_without_prior if
                                            d['marking'] == marking_process_net_dict]
    # example: [{'marking': {'p_1': 1}, 'outgoing_transitions':
    #                                        [{'unique_name': 't_A', 'frequency': 203, 'model_move_probability': 1.0}]}]

    # set of Transition objects, that are currently in the process net enabled
    enabled_transitions_for_marking = petri.semantics.enabled_transitions(process_net, marking_process_net)

    if len(model_move_probabilities_for_marking) == 1:
        # first, look if probability with applied prior has been already calculated
        for t in model_move_probabilities_for_marking[0]['outgoing_transitions']:
            if t["unique_name"] == wanted_transition and "model_move_probability_with_applied_prior" in t:
                return t["model_move_probability_with_applied_prior"]

        # current marking exists in model_move_probabilities
        frequency_of_outgoing_arcs_for_marking = 0
        for e_t in enabled_transitions_for_marking:
            # e_t is a Transition object
            transition_found = False

            for t in model_move_probabilities_for_marking[0]['outgoing_transitions']:
                # t is a dict, e.g. {'unique_name': 't_A', 'frequency': 204, 'model_move_probability': 1.0}
                if t['unique_name'] == e_t.name:
                    transition_found = True
                    t['frequency'] += constant_prior
                    frequency_of_outgoing_arcs_for_marking += t['frequency']

            if not transition_found:
                # transition given the marking was never executed in the training set (model_move_probabilities)
                new = {'unique_name': e_t.name, 'frequency': constant_prior,
                       'model_move_probability_with_applied_prior': None}
                model_move_probabilities_for_marking[0]['outgoing_transitions'].append(new)
                frequency_of_outgoing_arcs_for_marking += constant_prior

        res = None
        # calculate possibilities for model move given the marking based on new (added prior) frequencies
        for t in model_move_probabilities_for_marking[0]['outgoing_transitions']:
            t['model_move_probability_with_applied_prior'] = t['frequency'] / frequency_of_outgoing_arcs_for_marking
            print("probability calculation of: " + str(t['unique_name']))
            print("for marking: " + str(marking_process_net))
            print(str(t['frequency']) + " / " + str(frequency_of_outgoing_arcs_for_marking))
            print(str(t['model_move_probability_with_applied_prior']) + "\n")
            if t['unique_name'] == wanted_transition:
                res = t['model_move_probability_with_applied_prior']
        print(wanted_transition)
        print(marking_process_net_dict)
        print(str(res) + "\n")
        return res

    elif len(model_move_probabilities_for_marking) == 0:
        # current marking does not exist in model_move_probabilities
        frequency_of_outgoing_arcs_for_marking = 0
        unseen_marking = {'marking': marking_process_net_dict, 'outgoing_transitions': []}
        for e_t in enabled_transitions_for_marking:
            # e_t is a Transition object
            frequency_of_outgoing_arcs_for_marking += constant_prior
            unseen_marking['outgoing_transitions'].append(
                {'unique_name': e_t.name, 'frequency': constant_prior,
                 'model_move_probability_with_applied_prior': None})
        res = None
        # calculate possibilities for model move
        for t in unseen_marking['outgoing_transitions']:
            p = t['frequency'] / frequency_of_outgoing_arcs_for_marking
            t['model_move_probability_with_applied_prior'] = p
            if t['unique_name'] == wanted_transition:
                res = t['model_move_probability_with_applied_prior']
        model_move_probabilities_without_prior.append(unseen_marking)
        print("unseen marking")
        print(wanted_transition)
        print(marking_process_net_dict)
        print(res)
        print()
        return res
    else:
        # something went wrong, markings should be unique
        raise Exception('Multiple markings founded')


def calculate_model_move_prior():
    pass


def __get_log_move_probability(transition, log_move_probabilities):
    return log_move_probabilities[transition.label[0]]


def __compute_exact_heuristic(sync_net, incidence_matrix, marking):
    # TODO
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
