from pm4py.objects.log.importer.xes import factory as xes_importer
from pm4py.algo.conformance import alignments as alignments_module

SKIP = '>>'
STD_MODEL_LOG_MOVE_COST = 10000
STD_TAU_COST = 1
STD_SYNC_COST = 0


def calculate_log_move_probability(event_log, prior=None):
    event_frequencies_log = count_event_frequencies(event_log)
    event_frequencies = sum_up_prior_and_log_event_frequencies(event_frequencies_log, prior)
    total_number_event_names = len(event_frequencies)
    total_number_events = 0
    for event_name in event_frequencies:
        total_number_events += event_frequencies[event_name]
    log_move_probability = {**event_frequencies}
    for event_name in log_move_probability:
        log_move_probability[event_name] = (1 - event_frequencies[event_name] / total_number_events) / (
                total_number_event_names - 1)
    return log_move_probability


def count_event_frequencies(event_log):
    event_frequencies = {}
    for trace in event_log:
        for event in trace.get_event_list():
            event_name = event.get_dict()['concept:name']
            if event_name in event_frequencies:
                event_frequencies[event_name] += 1
            else:
                event_frequencies[event_name] = 1
    return event_frequencies


def sum_up_prior_and_log_event_frequencies(event_frequencies_log, event_frequencies_prior):
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

    for alignment in alignments:
        print("--------------------------")
        for move in alignment['alignment']:
            print("Marking before:  " + str(move['marking_before_transition']))
            print("Label:           " + str(move['label']))
            print("Name:            " + str(move['name']))
            print("Marking after:   " + str(move['marking_after_transition']))
            print()

