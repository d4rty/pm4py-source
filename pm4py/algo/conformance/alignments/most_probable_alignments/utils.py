from pm4py.algo.conformance import alignments as alignments_module
from pm4py.algo.conformance.alignments.utils import SKIP  # skip symbol, usual: '>>'


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

    # if parameters are none, the standard cost function will be applied ->
    #   cost 1000   log/model moves
    #   cost 1      tau moves
    #   cost 0      synchronous moves
    alignments = alignments_module.factory.apply(event_log, petri_net, initial_marking, final_marking, parameters=None)

    log_move_probabilities_given_marking = []

    for alignment in alignments:
        for move in alignment['alignment']:
            # only consider synchronous and model moves
            if (move['label'][0] == SKIP and move['label'][1] != SKIP) or (move['label'][0] == move['label'][1]):
                current_marking = {}
                for place in move['marking_before_transition']:
                    # example place name of the log net ('p_2','>>')
                    # example place name of the process net ('>>','p_2')
                    if place.name[0] == SKIP and place.name[1] != SKIP:
                        # only consider process net places
                        number_tokens = move['marking_before_transition'][place]
                        current_marking[place.name[1]] = number_tokens

                founded_dict = next(
                    (x for x in log_move_probabilities_given_marking if x['marking'] == current_marking), None)

                if founded_dict:
                    founded_transition = next(
                        (x for x in founded_dict['outgoing_transitions'] if
                         str(x['unique_name']) == str(move['name'][1])),
                        None)

                    if founded_transition:
                        founded_transition["frequency"] += 1
                    else:
                        new_transition = {"unique_name": str(move['name'][1]), "frequency": 1}
                        founded_dict["outgoing_transitions"].append(new_transition)
                else:
                    new_marking = {
                        "marking": current_marking,
                        "outgoing_transitions": [
                            {"unique_name": str(move['name'][1]), "frequency": 1, "model_move_probability": None}]
                    }
                    log_move_probabilities_given_marking.append(new_marking)

    for marking in log_move_probabilities_given_marking:
        total_number_fired_transitions_from_marking = 0
        for transition in marking['outgoing_transitions']:
            total_number_fired_transitions_from_marking += transition["frequency"]
        for transition in marking['outgoing_transitions']:
            transition["model_move_probability"] = transition["frequency"] / total_number_fired_transitions_from_marking
    return log_move_probabilities_given_marking


def print_most_probable_alignment(alignment):
    """
    Takes an alignment and prints the alignment to the console, e.g.:
     A  | B  | C  | D  |
    --------------------
     A  | B  | C  | >> |

    :param alignment: <class 'dict'>
    :return: Nothing
    """
    print("\nprobability: ", alignment['probability'], ' ~%.2f' % (alignment['probability'] * 100), '%')

    trace_steps = []
    model_steps = []
    max_label_length = 0
    for step in alignment['alignment']:
        trace_steps.append(" " + str(step['label'][0]) + " ")
        model_steps.append(" " + str(step['label'][1]) + " ")
        if len(step['label'][0]) > max_label_length:
            max_label_length = len(str(step['label'][0]))
        if len(str(step['label'][1])) > max_label_length:
            max_label_length = len(str(step['label'][1]))
    for i in range(len(trace_steps)):
        if len(str(trace_steps[i])) - 2 < max_label_length:
            step_length = len(str(trace_steps[i])) - 2
            spaces_to_add = max_label_length - step_length
            for j in range(spaces_to_add):
                if j % 2 == 0:
                    trace_steps[i] = trace_steps[i] + " "
                else:
                    trace_steps[i] = " " + trace_steps[i]
        print(trace_steps[i], end='|')
    divider = ""
    length_divider = len(trace_steps) * (max_label_length + 3)
    for i in range(length_divider):
        divider += "-"
    print('\n' + divider)
    for i in range(len(model_steps)):
        if len(model_steps[i]) - 2 < max_label_length:
            step_length = len(model_steps[i]) - 2
            spaces_to_add = max_label_length - step_length
            for j in range(spaces_to_add):
                if j % 2 == 0:
                    model_steps[i] = model_steps[i] + " "
                else:
                    model_steps[i] = " " + model_steps[i]

        print(model_steps[i], end='|')
    print('\n')
