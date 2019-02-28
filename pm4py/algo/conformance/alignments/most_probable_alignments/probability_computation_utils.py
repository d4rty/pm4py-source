import math

from pm4py.algo.conformance import alignments as alignments_module
from pm4py.algo.conformance.alignments.most_probable_alignments.miscellaneous_utils import \
    sync_product_net_place_belongs_to_process_net, is_model_move, is_log_move
from pm4py.algo.conformance.alignments.utils import SKIP  # skip symbol, usual: '>>'
from pm4py.algo.filtering.tracelog.variants import variants_filter
from pm4py.objects import petri
from pm4py.objects.petri.petrinet import Marking
from pm4py.statistics.traces.tracelog import case_statistics


def calculate_log_move_probability(event_log, prior=None):
    """
    Given an event_log and an optional prior the method calculates the log move probabilities for each unique event
    :param event_log: <class 'pm4py.objects.log.log.TraceLog'>
    :param prior: <class 'dict'>, e.g.: {'A': 1, 'B': 1, 'C': 1, 'D': 1, '*': 1}
    :return: <class 'dict'>, e.g.: {'A': 0.18795620437956204, 'B': 0.12439172749391728, ...}
    """
    event_frequencies_log = _count_event_frequencies(event_log)
    event_frequencies = _sum_up_prior_and_log_event_frequencies(event_frequencies_log, prior)
    total_number_event_names = len(event_frequencies)
    total_number_events = 0
    for event_name in event_frequencies:
        total_number_events += event_frequencies[event_name]
    log_move_probability = {**event_frequencies}
    for event_name in log_move_probability:
        log_move_probability[event_name] = (1 - event_frequencies[event_name] / total_number_events) / (
                total_number_event_names - 1)
    return log_move_probability


def _count_event_frequencies(event_log):
    """
    Takes a log and counts the occurrence of each unique event
    :param event_log: <class 'pm4py.objects.log.log.TraceLog'>
    :return: <class 'dict'>, e.g.: {'A': 203, 'B': 412, 'C': 162, 'D': 40}
    """
    event_frequencies = {}
    for trace in event_log:
        for event in trace.get_event_list():
            event_name = event.get_dict()['concept:name']
            if event_name in event_frequencies:
                event_frequencies[event_name] += 1
            else:
                event_frequencies[event_name] = 1
    return event_frequencies


def _sum_up_prior_and_log_event_frequencies(event_frequencies_log, event_frequencies_prior):
    """
    :param event_frequencies_log: <class 'dict'>, e.g.: {'A': 203, 'B': 412, 'C': 162, 'D': 40}
    :param event_frequencies_prior: <class 'dict'>, e.g.: {'A': 1, 'B': 1, 'C': 1, 'D': 1, '*': 1}
    :return: <class 'dict'>, e.g.: {'A': 204, 'B': 413, 'C': 163, 'D': 41, '*': 1}
    """
    res = {}
    if event_frequencies_prior:
        for event in event_frequencies_log:
            if event in event_frequencies_prior:
                res[event] = event_frequencies_prior[event] + event_frequencies_log[event]
            else:
                res[event] = event_frequencies_log[event]
        res = {**event_frequencies_prior, **res}
        return res
    else:
        return event_frequencies_log


def calculate_model_move_probabilities_without_prior(event_log, petri_net, initial_marking, final_marking):
    """
    Function calculates the model move probabilities for a petri net (with initial/final marking) and an event log
    :param event_log: <class 'pm4py.objects.log.log.TraceLog'>
    :param petri_net: <class 'pm4py.objects.petri.petrinet.PetriNet'>
    :param initial_marking: <class 'pm4py.objects.petri.petrinet.Marking'>
    :param final_marking: <class 'pm4py.objects.petri.petrinet.Marking'>
    :return: <class 'list'>, e.g.:
    [{'marking': {'p_1': 1}, 'outgoing_transitions':
                                            [{'unique_name': 't_A', 'frequency': 203, 'model_move_probability': 1.0}]
     },...]
    """

    variants_count = case_statistics.get_variant_statistics(event_log)

    all_opt_align_for_all_variants = []
    variants = variants_filter.get_variants(event_log)

    for unique_trace in variants:
        all_opt_align_trace = alignments_module.factory.apply(variants[unique_trace][0],
                                                              petri_net, initial_marking,
                                                              final_marking,
                                                              parameters=None,
                                                              find_all_opt_alignments=True)
        trace_frequency = next((variant["count"] for variant in variants_count if variant['variant'] == unique_trace),
                               None)
        all_opt_align_trace["trace"] = unique_trace
        all_opt_align_trace["trace_frequency"] = trace_frequency
        all_opt_align_for_all_variants.append(all_opt_align_trace)

    # if parameters are none, the standard cost function will be applied ->
    #   cost 1000   log/model moves
    #   cost 1      tau moves
    #   cost 0      synchronous moves

    model_move_probabilities_given_marking = []

    for trace in all_opt_align_for_all_variants:
        number_optimal_alignments = len(trace['alignments'])
        for alignment in trace['alignments']:
            for move in alignment:
                # only consider synchronous and model moves
                if (move['label'][0] == SKIP and move['label'][1] != SKIP) or (
                        move['label'][0] == move['label'][1]):
                    current_marking = {}
                    for place in move['marking_before_transition']:
                        # example place name of the log net ('p_2','>>')
                        # example place name of the process net ('>>','p_2')
                        if place.name[0] == SKIP and place.name[1] != SKIP:
                            # only consider process net places
                            number_tokens = move['marking_before_transition'][place]
                            current_marking[place.name[1]] = number_tokens

                    founded_dict = next(
                        (x for x in model_move_probabilities_given_marking if x['marking'] == current_marking),
                        None)

                    if founded_dict:
                        founded_transition = next(
                            (x for x in founded_dict['outgoing_transitions'] if
                             str(x['unique_name']) == str(move['name'][1])), None)
                        if founded_transition:
                            # increase number of already known transition given the marking
                            founded_transition["frequency"] += (1 / number_optimal_alignments) * trace[
                                "trace_frequency"]
                        else:
                            # add new transition to already known marking
                            new_transition = {"unique_name": str(move['name'][1]),
                                              "frequency": (1 / number_optimal_alignments) * trace["trace_frequency"]}
                            founded_dict["outgoing_transitions"].append(new_transition)
                    else:
                        new_marking = {
                            "marking": current_marking,
                            "outgoing_transitions": [
                                {"unique_name": str(move['name'][1]),
                                 "frequency": (1 / number_optimal_alignments) * trace["trace_frequency"],
                                 "model_move_probability": None}]
                        }
                        model_move_probabilities_given_marking.append(new_marking)

    for marking in model_move_probabilities_given_marking:
        total_number_fired_transitions_from_marking = 0
        for transition in marking['outgoing_transitions']:
            total_number_fired_transitions_from_marking += transition["frequency"]
        for transition in marking['outgoing_transitions']:
            transition["model_move_probability"] = transition[
                                                       "frequency"] / total_number_fired_transitions_from_marking

    return model_move_probabilities_given_marking


def get_move_probability(transition, marking, log_move_probabilities, model_move_probabilities, process_net):
    """

    :param transition:
    :param marking:
    :param log_move_probabilities:
    :param model_move_probabilities:
    :param process_net:
    :return:
    """
    if is_model_move(transition, SKIP):
        return get_model_move_probability(transition, marking, model_move_probabilities, process_net)
    elif is_log_move(transition, SKIP):
        return get_log_move_probability(transition, log_move_probabilities)
    else:
        # synchronous move
        cost = max(get_log_move_probability(transition, log_move_probabilities),
                   get_model_move_probability(transition, marking, model_move_probabilities, process_net))
        return cost


def get_log_move_probability(transition, log_move_probabilities):
    return log_move_probabilities[transition.label[0]]


def get_model_move_probability(transition, marking_sync_product_net, model_move_probabilities_without_prior,
                               process_net, constant_prior=0):
    """
    :param transition: Transition object
    :param marking_sync_product_net: Marking object
    :param model_move_probabilities_without_prior: list of dicts (key: marking, value: outgoing transitions and
                                                                                                    their frequencies)
    :param process_net: PetriNet object
    :param constant_prior: Prior that is added to all transitions
    :return: Double - probability of the transition given the marking
    """
    wanted_transition = transition.name[1]

    # retrieve the marking from the process net part of the synchronous product net
    marking_process_net = Marking()
    marking_process_net_as_dict = {}  # dict of marking - key: place name, value: number tokens

    for place in marking_sync_product_net:
        if sync_product_net_place_belongs_to_process_net(place):
            number_tokens = marking_sync_product_net[place]
            process_net_place_name = place.name[1]
            # create marking object of process net
            for p in process_net.places:
                if p.name == process_net_place_name:
                    marking_process_net[p] = number_tokens
            marking_process_net_as_dict[process_net_place_name] = number_tokens

    model_move_probabilities_for_marking = [d for d in model_move_probabilities_without_prior if
                                            d['marking'] == marking_process_net_as_dict]
    # example of model_move_probabilities_for_markingib: [{'marking': {'p_1': 1}, 'outgoing_transitions':
    #                                        [{'unique_name': 't_A', 'frequency': 203, 'model_move_probability': 1.0}]}]

    # set of Transition objects, that are currently in the process net enabled
    enabled_transitions_for_marking = petri.semantics.enabled_transitions(process_net, marking_process_net)

    if len(model_move_probabilities_for_marking) == 1:
        # case 1: current marking of process net exists in model_move_probabilities

        # first, look if probability with applied prior has been already calculated
        for t in model_move_probabilities_for_marking[0]['outgoing_transitions']:
            if t["unique_name"] == wanted_transition and "model_move_probability_with_applied_prior" in t:
                return t["model_move_probability_with_applied_prior"]

        frequency_of_outgoing_arcs_for_marking = 0

        # holds transitions that have not been executed before (given the marking) in the training set
        new_transitions = []
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
                new_transitions.append(new)
                frequency_of_outgoing_arcs_for_marking += constant_prior

        model_move_probabilities_for_marking[0]['outgoing_transitions'].extend(new_transitions)
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
        return res

    elif len(model_move_probabilities_for_marking) == 0:
        # case 2: current marking of process net does not exist in model_move_probabilities
        frequency_of_outgoing_arcs_for_marking = 0
        unseen_marking = {'marking': marking_process_net_as_dict, 'outgoing_transitions': []}

        for e_t in enabled_transitions_for_marking:
            # e_t is a Transition object
            frequency_of_outgoing_arcs_for_marking += constant_prior
            unseen_marking['outgoing_transitions'].append(
                {'unique_name': e_t.name, 'frequency': constant_prior,
                 'model_move_probability_with_applied_prior': None})
        res = None
        # calculate possibilities for model move
        for t in unseen_marking['outgoing_transitions']:
            if frequency_of_outgoing_arcs_for_marking > 0:
                p = t['frequency'] / frequency_of_outgoing_arcs_for_marking
            else:
                p = 0
            t['model_move_probability_with_applied_prior'] = p
            if t['unique_name'] == wanted_transition:
                res = t['model_move_probability_with_applied_prior']
        model_move_probabilities_without_prior.append(unseen_marking)

        print("unseen marking")
        print(wanted_transition)
        print(marking_process_net_as_dict)
        print(res)
        print()
        return res
    else:
        # something went wrong, markings should be unique
        raise Exception('Multiple markings founded')


def apply_log_transformation(probability):
    if probability is None:
        raise ValueError('Probability is of type NoneType')
    # TODO optimize --> apply log transformation not on the fly but apply it if possible before
    # --> try to avoid duplicate calls
    if not (0 <= probability <= 1):
        raise ValueError('probability needs to be in the range [0,1]')
    if probability > 0:
        return - math.log(probability)
    elif probability == 0:
        # -log(0) ~ math.inf
        return math.inf
    else:
        raise ValueError('Log Transformation not applicable to numbers below 0')


def get_move_cost(transition, marking, log_move_probabilities, model_move_probabilities, process_net):
    """

    :param transition:
    :param marking:
    :param log_move_probabilities:
    :param model_move_probabilities:
    :param process_net:
    :return:
    """
    if is_model_move(transition, SKIP):
        return apply_log_transformation(
            get_model_move_probability(transition, marking, model_move_probabilities, process_net))
    elif is_log_move(transition, SKIP):
        return apply_log_transformation(get_log_move_probability(transition, log_move_probabilities))
    else:
        # synchronous move
        cost = min(
            apply_log_transformation(get_log_move_probability(transition, log_move_probabilities)),
            apply_log_transformation(
                get_model_move_probability(transition, marking, model_move_probabilities, process_net)))
        return cost
