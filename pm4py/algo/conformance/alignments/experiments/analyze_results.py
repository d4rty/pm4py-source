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
    plot_time_per_algorithm, generate_simple_bar_plot


def start():
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
        'Incremental\nA* with\nheuristic',
        'online\nconformance\nno window size',
        'online\nconformance\nwindow size: 1',
        'online\nconformance\nwindow size: 2',
        'online\nconformance\nwindow size: 3',
    ]

    path_to_files = os.path.join("C:\\", "Users", "Daniel", "Desktop", "master_thesis", "experiments",
                                 "results", 'bpi_ch_19')
    files = os.listdir(path_to_files)
    result_files = [
        'RESULTS_petri_net_1.pnml_2019-04-27.pickle',
        'RESULTS_petri_net_2.pnml_2019-04-27.pickle',
        'RESULTS_petri_net_3.pnml_2019-04-27.pickle',
        'RESULTS_petri_net_4.pnml_2019-04-27.pickle',
        'RESULTS_petri_net_5.pnml_2019-04-27.pickle'
    ]

    for res_file in result_files:
        print(res_file)
        pickle_path = os.path.join(path_to_files, res_file)
        with open(pickle_path, 'rb') as handle:
            results = pickle.load(handle)

        plot_time_bar_chart(algo_result_keys, results, path_to_files, res_file, description_algos)
        bar_plot_miscellaneous(algo_result_keys, results, path_to_files, res_file, description_algos, 'traversed_arcs')
        bar_plot_miscellaneous(algo_result_keys, results, path_to_files, res_file, description_algos, 'visited_states')
        bar_plot_miscellaneous(algo_result_keys, results, path_to_files, res_file, description_algos, 'queued_states')


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
    filename = os.path.join(path_to_files, attribute + '_' + res_file)
    generate_simple_bar_plot(tuple(res), tuple(description_algos), filename, svg=False)

    print("plot " + attribute + " for a subset of the algorithm variants")
    # time plot for subset of algorithm variants
    subset_variants = [2, 4, 5, 6, 7, 8]
    res = [res[i - 1] for i in
           subset_variants]

    description_algos_subset = [description_algos[i - 1] for i in subset_variants]

    filename = os.path.join(path_to_files, attribute + '_subset_' + res_file)
    generate_simple_bar_plot(tuple(res), tuple(description_algos_subset), filename, svg=False)


if __name__ == '__main__':
    start()
