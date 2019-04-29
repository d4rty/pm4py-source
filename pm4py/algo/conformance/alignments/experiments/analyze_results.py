import datetime
import os
import pickle
from datetime import date
from pm4py.objects import petri
from pm4py.algo.conformance.alignments.incremental_a_star.incremental_a_star import \
    apply as incremental_a_star_apply
import pandas as pd
from random import *
from pm4py.algo.conformance.alignments.versions.state_equation_a_star import apply as state_equation_a_star_apply
from pm4py.algo.conformance.alignments.incremental_a_star.incremental_prefix_alignments_bas import \
    apply as incremental_prefix_alignments_apply
from pm4py.objects.log.log import Trace
from pm4py.algo.conformance.alignments.experiments.plot import plot_length_distribution, \
    plot_time_per_algorithm, generate_simple_bar_plot, plot_search_space_size_depending_on_prefix_length


def start_bpi_2019():
    algo_result_keys = [
        'a_star_from_scratch_without_heuristic',
        'a_star_from_scratch_with_heuristic',
        'incremental_a_star_without_heuristic',
        'incremental_a_star_with_heuristic',
        'online_conformance_window_0',
        'online_conformance_window_1',
        'online_conformance_window_2',
        'online_conformance_window_3',
    ]

    description_algos = [
        'A* no\nheuristic',
        'A* with\nheuristic',
        'incremental\nA* no\n heuristic',
        'incremental\nA* with\nheuristic',
        'online\nconformance\nno window size',
        'online\nconformance\nwindow size: 1',
        'online\nconformance\nwindow size: 2',
        'online\nconformance\nwindow size: 3',
    ]

    algo_result_keys_to_description_string = {
        'a_star_from_scratch_without_heuristic': 'A* no\nheuristic',
        'a_star_from_scratch_with_heuristic': 'A* with\nheuristic',
        'incremental_a_star_without_heuristic': 'incremental\nA* no\n heuristic',
        'incremental_a_star_with_heuristic': 'incremental\nA* with\nheuristic',
        'online_conformance_window_0': 'online\nconformance\nno window size',
        'online_conformance_window_1': 'online\nconformance\nwindow size: 1',
        'online_conformance_window_2': 'online\nconformance\nwindow size: 2',
        'online_conformance_window_3': 'online\nconformance\nwindow size: 3',
    }

    path_to_files = os.path.join("C:\\", "Users", "Daniel", "Desktop", "master_thesis", "experiments",
                                 "results", 'bpi_ch_19')
    result_files = [
        'RESULTS_petri_net_1.pnml_2019-04-27.pickle',
        'RESULTS_petri_net_2.pnml_2019-04-27.pickle',
        'RESULTS_petri_net_3.pnml_2019-04-27.pickle',
        'RESULTS_petri_net_4.pnml_2019-04-27.pickle',
        'RESULTS_petri_net_5.pnml_2019-04-27.pickle'
    ]

    measured_attributes = ['traversed_arcs', 'visited_states', 'queued_states', 'total_computation_time']

    for res_file in result_files:
        print(res_file)
        pickle_path = os.path.join(path_to_files, res_file)
        with open(pickle_path, 'rb') as handle:
            results = pickle.load(handle)

        # just to check if costs are the same #########################################################################
        # for trace in results:
        #     cost = None
        #     for key in algo_result_keys[:5]:
        #         if cost:
        #             if not int(trace[key]['cost']/1000) == cost:
        #                 print("Error")
        #         else:
        #             cost = int(trace[key]['cost']/1000)

        plot_time_bar_chart(algo_result_keys, results, path_to_files, res_file, description_algos)
        bar_plot_miscellaneous(algo_result_keys, results, path_to_files, res_file, description_algos, 'traversed_arcs')
        bar_plot_miscellaneous(algo_result_keys, results, path_to_files, res_file, description_algos, 'visited_states')
        bar_plot_miscellaneous(algo_result_keys, results, path_to_files, res_file, description_algos, 'queued_states')

        # remove line breaks
        for k in algo_result_keys_to_description_string:
            algo_result_keys_to_description_string[k] = algo_result_keys_to_description_string[k].replace("\n", " ")

        for attr in measured_attributes:
            store_plot_path = os.path.join(path_to_files, res_file + '_' + attr + '_given_prefix_length')
            res = get_attribute_per_prefix_length(algo_result_keys, results, description_algos, attr)
            plot_search_space_size_depending_on_prefix_length(algo_result_keys, algo_result_keys_to_description_string,
                                                              res, attr, path_to_store=store_plot_path, svg=False)
            del res['a_star_from_scratch_without_heuristic']
            del res['incremental_a_star_without_heuristic']
            store_plot_path = os.path.join(path_to_files, res_file + '_' + attr + '_given_prefix_length_SUBSET')
            plot_search_space_size_depending_on_prefix_length(algo_result_keys, algo_result_keys_to_description_string,
                                                              res, attr, path_to_store=store_plot_path, svg=False)


def plot_time_bar_chart(algo_result_keys, results, path_to_files, res_file, description_algos):
    a_star_computation_time_without_heuristic = [0] * len(algo_result_keys)
    computation_time_heuristic = [0] * len(algo_result_keys)
    for trace in results:
        for i, algo_variant in enumerate(algo_result_keys):
            time_without_heuristic = trace[algo_variant]['total_computation_time'] - \
                                     trace[algo_variant]['heuristic_computation_time']
            a_star_computation_time_without_heuristic[i] += time_without_heuristic
            computation_time_heuristic[i] += trace[algo_variant]['heuristic_computation_time']

    a_star_computation_time_without_heuristic = [0.01 * v for v in a_star_computation_time_without_heuristic]
    computation_time_heuristic = [0.01 * v for v in computation_time_heuristic]

    print("plot time for all algorithm variants")
    # time plot for all algorithm variants
    filename = os.path.join(path_to_files, 'computation_time_all_' + res_file)
    plot_time_per_algorithm(tuple(a_star_computation_time_without_heuristic), tuple(computation_time_heuristic),
                            tuple(description_algos), filename, svg=False)

    print("plot time for subset of algorithm variants")
    # time plot for subset of algorithm variants
    subset_variants = [2, 4, 5, 6, 7, 8]
    a_star_computation_time_without_heuristic = [a_star_computation_time_without_heuristic[i - 1] for i in
                                                 subset_variants]
    computation_time_heuristic = [computation_time_heuristic[i - 1] for i in subset_variants]
    description_algos_subset = [description_algos[i - 1] for i in subset_variants]

    filename = os.path.join(path_to_files, 'computation_time_subset_' + res_file)
    plot_time_per_algorithm(tuple(a_star_computation_time_without_heuristic),
                            tuple(computation_time_heuristic),
                            tuple(description_algos_subset), filename, svg=False)


def bar_plot_miscellaneous(algo_result_keys, results, path_to_files, res_file, description_algos, attribute):
    res = [0] * len(algo_result_keys)
    for trace in results:
        for i, algo_variant in enumerate(algo_result_keys):
            res[i] += trace[algo_variant][attribute]

    res = [0.01 * v for v in res]  # compute average

    print("plot " + attribute + " for all algorithm variants")
    # time plot for all algorithm variants
    filename = os.path.join(path_to_files, res_file + '_' + attribute)
    generate_simple_bar_plot(tuple(res), tuple(description_algos), filename, svg=False)

    print("plot " + attribute + " for a subset of the algorithm variants")
    # time plot for subset of algorithm variants
    subset_variants = [2, 4, 5, 6, 7, 8]
    res = [res[i - 1] for i in
           subset_variants]

    description_algos_subset = [description_algos[i - 1] for i in subset_variants]

    filename = os.path.join(path_to_files, res_file + '_' + attribute + '_subset')
    generate_simple_bar_plot(tuple(res), tuple(description_algos_subset), filename, svg=False)


def get_attribute_per_prefix_length(algo_result_keys, results, description_algos, attribute):
    res = {}
    for j, algorithm in enumerate(algo_result_keys):
        number_prefix_lengths = {}  # how often was a certain prefix lengths in intermediate results
        res[algorithm] = []

        for trace in results:

            for i, intermediate_result in enumerate(trace[algorithm]["intermediate_results"]):

                if i + 1 in number_prefix_lengths:
                    number_prefix_lengths[i + 1] += 1
                else:
                    number_prefix_lengths[i + 1] = 1

                if len(res[algorithm]) <= i:
                    res[algorithm].append(intermediate_result[attribute])
                else:
                    res[algorithm][i] += intermediate_result[attribute]

        # calculate cumulated numbers per prefix length
        for i in range(1, len(res[algorithm])):
            res[algorithm][i] = res[algorithm][i] + res[algorithm][i - 1]

        # calculate average
        # for i in range(len(res[algorithm])):
        #     res[algorithm][i] = (1 / number_prefix_lengths[i + 1]) * res[algorithm][i]
    return res


if __name__ == '__main__':
    start_bpi_2019()
