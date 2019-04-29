import datetime
import os
import pickle
from datetime import date

from pm4py.algo.conformance.alignments.experiments.get_search_space_size_info import plot_search_space_size
from pm4py.objects.log.importer.xes import factory as xes_importer
from pm4py.objects import petri
from pm4py.algo.conformance.alignments.incremental_a_star.incremental_a_star import \
    apply as incremental_a_star_apply
import pandas as pd
from random import *
from pm4py.algo.conformance.alignments.versions.state_equation_a_star import apply as state_equation_a_star_apply
from pm4py.algo.conformance.alignments.incremental_a_star.incremental_prefix_alignments_bas import \
    apply as incremental_prefix_alignments_apply
from pm4py.objects.log.log import Trace
from pm4py.algo.conformance.alignments.experiments.plot import plot_length_distribution


def calculate_prefix_alignments(petri_net_filename, log, path_to_files, trace_indices):
    results = []

    pnml_file_path = os.path.join(path_to_files, petri_net_filename)
    net, im, fm = petri.importer.pnml.import_net(
        pnml_file_path)
    # gviz = pn_vis_factory.apply(net, marking, fmarking,
    #                             parameters={"format": "svg", 'debug': False})
    # pn_vis_factory.view(gviz)

    iteration = 1
    for i in trace_indices:
        print("Trace index: %i | Trace %i of %i | %s\n" % (i, iteration, len(trace_indices), petri_net_filename))
        print(log[i])

        print("calculate_prefix_alignments_dijkstra_from_scratch")
        res0 = calculate_prefix_alignments_from_scratch(log[i], net, im, fm, dijkstra=True)

        print("calculate_prefix_alignments_from_scratch_with_heuristic")
        res1 = calculate_prefix_alignments_from_scratch(log[i], net, im, fm, dijkstra=False)

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

        print("calculate_prefix_alignment_online_conformance_by_bas - window size 3")
        res7 = calculate_prefix_alignment_online_conformance_by_bas(log[i], net, im, fm, window_size=3)

        results.append({'trace': log[i],
                        'trace_index': i,
                        'a_star_from_scratch_without_heuristic': res0,
                        'a_star_from_scratch_with_heuristic': res1,
                        'incremental_a_star_without_heuristic': res2,
                        'incremental_a_star_with_heuristic': res3,
                        'online_conformance_window_0': res4,
                        'online_conformance_window_1': res5,
                        'online_conformance_window_2': res6,
                        'online_conformance_window_3': res7
                        })
        iteration += 1
    results_path = os.path.join(path_to_files, "RESULTS_" + petri_net_filename + '_' + str(date.today()) + ".pickle")
    with open(results_path, 'wb') as handle:
        pickle.dump(results, handle)
    # with open(results_path, 'rb') as handle:
    #     b = pickle.load(handle)
    #     print(b)


def calculate_prefix_alignments_from_scratch(trace, petri_net, initial_marking, final_marking, dijkstra):
    '''
    This method calculates for a trace prefix alignments by starting a* WITHOUT a heurisitc (i.e. dijkstra) every time
    from scratch, e.g, for <e_1, e2, .. e_n>, this methods calculates first an alignment for <e_1>,
    afterwards <e_1,e_2>, ... and so on
    :return:
    '''

    visited_states_total = 0
    queued_states_total = 0
    traversed_arcs_total = 0
    total_computation_time_total = 0
    heuristic_computation_time_total = 0
    number_solved_lps_total = 0

    intermediate_results = []
    incremental_trace = Trace()
    for event in trace:
        incremental_trace.append(event)
        res = state_equation_a_star_apply(incremental_trace, petri_net, initial_marking, final_marking,
                                          dijkstra=dijkstra)

        intermediate_result = {'trace_length': len(incremental_trace),
                               'alignment': res['alignment'],
                               'cost': res['cost'],
                               'visited_states': res['visited_states'],
                               'queued_states': res['queued_states'],
                               'traversed_arcs': res['traversed_arcs'],
                               'total_computation_time': res['total_computation_time'],
                               'heuristic_computation_time': res['heuristic_computation_time'],
                               'number_solved_lps': res['number_solved_lps']}

        visited_states_total += res['visited_states']
        queued_states_total += res['queued_states']
        traversed_arcs_total += res['traversed_arcs']
        total_computation_time_total += res['total_computation_time']
        heuristic_computation_time_total += res['heuristic_computation_time']
        number_solved_lps_total += res['number_solved_lps']

        intermediate_results.append(intermediate_result)
    res['intermediate_results'] = intermediate_results
    res['visited_states'] = visited_states_total
    res['queued_states'] = queued_states_total
    res['traversed_arcs'] = traversed_arcs_total
    res['total_computation_time'] = total_computation_time_total
    res['heuristic_computation_time'] = heuristic_computation_time_total
    res['number_solved_lps'] = number_solved_lps_total
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


def get_x_traces_with_min_length_from_log(log, number_traces, min_trace_length, max_length_trace):
    '''
    :param log:
    :return: array (of length number_traces) consisting indices that correspond to traces in the log that have at least
    min_trace_length events
    '''
    res = []
    number_traces_log = len(log)
    for i in range(number_traces):
        found_trace = False
        while not found_trace:
            trace_index = sample(range(number_traces_log), 1)[0]
            if min_trace_length <= len(log[trace_index]) <= max_length_trace and trace_index not in res:
                res.append(trace_index)
                found_trace = True
        print(i)
    return res


def execute_experiments_for_bpi_ch_2019():
    path_to_files = os.path.join("C:\\", "Users", "Daniel", "Desktop", "master_thesis", "experiments",
                                 "data", "bpi_challenge_2019")
    pickle_path = os.path.join(path_to_files, 'log.pickle')
    with open(pickle_path, 'rb') as handle:
        log = pickle.load(handle)

    # print(get_x_traces_with_min_length_from_log(log, 100, 1, 10))
    # see result below

    # indices of 100 traces that have a min length of 10
    trace_indices = [156987, 223136, 185014, 83947, 30839, 208130, 141950, 224143, 84577, 199169, 100002, 231178,
                     125148, 87677, 159548, 53192, 60566, 133799, 30757, 179923, 128658, 57995, 100238, 217181, 53090,
                     17284, 106604, 108682, 93742, 104953, 67881, 217093, 57514, 52539, 78006, 70889, 183598, 51575,
                     239034, 163527, 81664, 157315, 171228, 37554, 198622, 137706, 183177, 187943, 135636, 134117,
                     120926, 78350, 77484, 168340, 182852, 208529, 132746, 155675, 17121, 76351, 101531, 72736, 75676,
                     174786, 174541, 193458, 216456, 75885, 27423, 182877, 248427, 122730, 43213, 203694, 141650,
                     124713, 243940, 116283, 221914, 246544, 63639, 223875, 240937, 84138, 97163, 55977, 23651, 85792,
                     87578, 126468, 151616, 146372, 219597, 222486, 161946, 108353, 246970, 232496, 221165, 225121]

    # calculate_prefix_alignments("petri_net_1.pnml", log, path_to_files, trace_indices)
    # calculate_prefix_alignments("petri_net_2.pnml", log, path_to_files, trace_indices)
    # calculate_prefix_alignments("petri_net_3.pnml", log, path_to_files, trace_indices)
    # calculate_prefix_alignments("petri_net_4.pnml", log, path_to_files, trace_indices)
    # calculate_prefix_alignments("petri_net_5.pnml", log, path_to_files, trace_indices)

    # trace length distribution ########################################################################################
    # lengths = []
    # for i in trace_indices:
    #     lengths.append(len(log[i]))
    # plot_length_distribution(lengths, path_to_files)

    # search space size per prefix length ##############################################################################
    # petri_net_file_names = ["petri_net_1.pnml", "petri_net_2.pnml", "petri_net_3.pnml", "petri_net_4.pnml",
    #                         "petri_net_5.pnml"]
    # net_paths = [os.path.join(path_to_files, net) for net in petri_net_file_names]
    # plot_search_space_size(log, trace_indices, net_paths)


def execute_experiments_for_bpi_ch_2018():
    path_to_files = os.path.join("C:\\", "Users", "Daniel", "Desktop", "master_thesis", "experiments",
                                 "data", "bpi_challenge_2018")
    pickle_path = os.path.join(path_to_files, "log.pickle")

    # file_path = os.path.join(path_to_files, "BPI Challenge 2018.xes")
    # log = xes_importer.import_log(file_path)
    # with open(pickle_path, 'wb') as handle:
    #     pickle.dump(log, handle)

    pickle_path = os.path.join(path_to_files, 'log.pickle')
    with open(pickle_path, 'rb') as handle:
        log = pickle.load(handle)

    # log has only traces of length between 24 and 2973
    # print("start find 100 traces")
    # print(get_x_traces_with_min_length_from_log(log, 100, 24, 40))
    # see result below

    # indices of 100 traces that have a min length of 10
    trace_indices = [27909, 42420, 24453, 32291, 29211, 25260, 38148, 27479, 7054, 27847, 30589, 20880, 23365, 25161,
                     23005, 15787, 28539, 27052, 22954, 16573, 20099, 25897, 19182, 16115, 18730, 22390, 30036, 18155,
                     39391, 43390, 25572, 23872, 14757, 15745, 10200, 16240, 42863, 15422, 36625, 37472, 24537, 26687,
                     28163, 27094, 36804, 20035, 22184, 33856, 14416, 28623, 16407, 16486, 22860, 26027, 18287, 22410,
                     24929, 20987, 35663, 19309, 32769, 17993, 18705, 34057, 16645, 23740, 19396, 17223, 36506, 23049,
                     22581, 37107, 23199, 26068, 13340, 29224, 22082, 19488, 37813, 28704, 16550, 25944, 17738, 42180,
                     24396, 16374, 19975, 25860, 37457, 25044, 25318, 25249, 25444, 28317, 23170, 19211, 28044, 13496,
                     27627, 21050]

    # calculate_prefix_alignments("petri_net_1.pnml", log, path_to_files, trace_indices)
    # calculate_prefix_alignments("petri_net_2.pnml", log, path_to_files, trace_indices)
    # calculate_prefix_alignments("petri_net_3.pnml", log, path_to_files, trace_indices)
    # calculate_prefix_alignments("petri_net_4.pnml", log, path_to_files, trace_indices)
    # calculate_prefix_alignments("petri_net_5.pnml", log, path_to_files, trace_indices)

    lengths = []
    for i in trace_indices:
        lengths.append(len(log[i]))
    plot_length_distribution(lengths, path_to_files)


if __name__ == '__main__':
    execute_experiments_for_bpi_ch_2019()
