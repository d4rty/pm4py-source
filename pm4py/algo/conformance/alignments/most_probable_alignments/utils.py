from pm4py.algo.conformance import alignments as alignments_module
from pm4py.algo.conformance.alignments.utils import SKIP  # skip symbol, usual: '>>'


def calculate_log_move_probability(event_log, prior=None):
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


def calculate_model_move_probability(event_log, petri_net, initial_marking, final_marking):
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

                # log_move_probabilities_given_marking is a list of dict like below:
                # {
                #     "marking": current_marking,
                #     "outgoing_transitions": [
                #         {"unique_name": move['name'][1], "frequency": xx, "model_move_probability": xx}
                #     ]
                # }

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
