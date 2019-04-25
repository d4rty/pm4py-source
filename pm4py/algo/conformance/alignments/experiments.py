import datetime
import os
import pickle

from pm4py.objects import petri
from pm4py.algo.conformance.alignments.incremental_a_star.incremental_a_star import \
    apply as incremental_a_star_apply
import pandas as pd
from random import *
from pm4py.algo.conformance.alignments.versions.state_equation_a_star import apply as state_equation_a_star_apply
from pm4py.algo.conformance.alignments.incremental_a_star.incremental_prefix_alignments_bas import \
    apply as incremental_prefix_alignments_apply


def calculate_prefix_alignments(petri_net_filename, log, path_to_files, trace_indices):
    results = {'alignment': [], 'trace': [], 'trace_index': [], 'computation_duration_total': [],
               'computation_duration_lp_solving': [], 'size_closed_set': []}

    pnml_file_path = os.path.join(path_to_files, petri_net_filename)
    net, im, fm = petri.importer.pnml.import_net(
        pnml_file_path)
    # gviz = pn_vis_factory.apply(net, marking, fmarking,
    #                             parameters={"format": "svg", 'debug': False})
    # pn_vis_factory.view(gviz)

    iteration = 1
    for i in trace_indices:
        print("Trace index: %i | Trace %i of %i | %s\n" % (i, iteration, len(trace_indices), petri_net_filename))
        print("calculate_prefix_alignments_dijkstra_from_scratch")
        res0 = calculate_prefix_alignments_dijkstra_from_scratch(log[i], net, im, fm)
        print("calculate_prefix_alignments_from_scratch_with_heuristic")
        res1 = calculate_prefix_alignments_from_scratch_with_heuristic(log[i], net, im, fm)
        print("calculate_prefix_alignment_modified_a_star_dijkstra")
        res2 = calculate_prefix_alignment_modified_a_star_dijkstra(log[i], net, im, fm)
        print("calculate_prefix_alignment_modified_a_star_with_heuristic")
        res3 = calculate_prefix_alignment_modified_a_star_with_heuristic(log[i], net, im, fm)
        print("calculate_prefix_alignment_online_conformance_by_bas")
        res4 = calculate_prefix_alignment_online_conformance_by_bas(log[i], net, im, fm)
        print("calculate_prefix_alignment_online_conformance_by_bas - window size 1")
        res5 = calculate_prefix_alignment_online_conformance_by_bas(log[i], net, im, fm, window_size=1)
        print("calculate_prefix_alignment_online_conformance_by_bas - window size 2")
        res6 = calculate_prefix_alignment_online_conformance_by_bas(log[i], net, im, fm, window_size=2)
        iteration += 1


def calculate_prefix_alignments_dijkstra_from_scratch(trace, petri_net, initial_marking, final_marking):
    '''
    This method calculates for a trace prefix alignments by starting a* WITHOUT a heurisitc (i.e. dijkstra) every time
    from scratch, e.g, for <e_1, e2, .. e_n>, this methods calculates first an alignment for <e_1>,
    afterwards <e_1,e_2>, ... and so on
    :return:
    '''
    res = state_equation_a_star_apply(trace, petri_net, initial_marking, final_marking, dijkstra=True)
    return res


def calculate_prefix_alignments_from_scratch_with_heuristic(trace, petri_net, initial_marking, final_marking):
    '''
    This method calculates for a trace prefix alignments by starting a* WITH a heurisitc (i.e. state equation) every
    time from scratch, e.g, for <e_1, e2, .. e_n>, this methods calculates first an alignment for <e_1>,
    afterwards <e_1,e_2>, ... and so on
    :return:
    '''

    res = state_equation_a_star_apply(trace, petri_net, initial_marking, final_marking)
    return res


def calculate_prefix_alignment_modified_a_star_dijkstra(trace, petri_net, initial_marking, final_marking):
    '''
    This method calculates for a trace prefix alignments by starting a* WITHOUT a heurisitc (i.e. dijkstra) and keeps
    open and closed set in memory
    :return:
    '''
    res = incremental_a_star_apply(trace, petri_net, initial_marking, final_marking, dijkstra=True)
    return res


def calculate_prefix_alignment_modified_a_star_with_heuristic(trace, petri_net, initial_marking, final_marking):
    '''
    This method calculates for a trace prefix alignments by starting a* WITH a heurisitc (i.e. modified state equation)
    and keeps open and closed set in memory
    :return:
    '''
    res = incremental_a_star_apply(trace, petri_net, initial_marking, final_marking)
    return res


def calculate_prefix_alignment_online_conformance_by_bas(trace, petri_net, initial_marking, final_marking,
                                                         window_size=0):
    '''
    This methods uses Bas' method WITH optimality guarantees (i.e. no partial reverting)
    :return:
    '''
    res = incremental_prefix_alignments_apply(trace, petri_net, initial_marking, final_marking, window_size=window_size)
    return res


# def generate_petri_net(log, number_samples=None):
#     for i in range(10):
#         print("start sampling")
#         if number_samples:
#             sampled_log = sampling.sample(log, n=number_samples)
#         else:
#             sampled_log = log
#         print("start mining")
#         net, initial_marking, final_marking = inductive_miner.apply(sampled_log)
#         gviz = pn_vis_factory.apply(net, initial_marking, final_marking,
#                                     parameters={"format": "svg", 'debug': False})
#         pn_vis_factory.view(gviz)
#     return net, initial_marking, final_marking


def get_x_traces_with_min_length_from_log(log, number_traces, min_trace_length):
    '''
    :param log:
    :return: array (of length number_traces) consisting indices that correspond to traces in the log that have at least
    min_trace_length events
    '''
    res = []
    number_traces_log = len(log)
    print(number_traces_log)
    for i in range(number_traces):
        found_trace = False
        while not found_trace:
            trace_index = sample(range(number_traces_log), 1)[0]
            if len(log[trace_index]) >= min_trace_length and trace_index not in res:
                res.append(trace_index)
                found_trace = True
    return res


def execute_experiments_for_bpi_ch_2019():
    path_to_files = os.path.join("C:\\", "Users", "Daniel", "Desktop", "master_thesis", "experiments",
                                 "data", "bpi_challenge_2019")
    pickle_path = os.path.join(path_to_files, 'log.pickle')
    with open(pickle_path, 'rb') as handle:
        log = pickle.load(handle)

    # print(get_x_traces_with_min_length_from_log(log, 100, 10))
    # see result below

    # indices of 100 traces that have a min length of 10
    trace_indices = [9863]

    calculate_prefix_alignments("petri_net_2.pnml", log, path_to_files, trace_indices)


if __name__ == '__main__':
    execute_experiments_for_bpi_ch_2019()
