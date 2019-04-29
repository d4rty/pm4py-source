from pm4py.objects import petri
import os
from pm4py.objects.log.importer.xes import factory as xes_importer
from pm4py.algo.conformance.alignments.utils import SKIP, print_alignment
from pm4py.objects.petri.reachability_graph import construct_reachability_graph
from pm4py.objects.log.log import Trace


def plot_search_space_size(log, trace_indices, net_paths):
    for net in net_paths:
        net, i_marking, f_marking = petri.importer.pnml.import_net(net)
        for i in trace_indices:
            print(len(log[i]))
            trace = Trace()
            for j in range(len(log[i])):
                trace.append(log[i][j])
                print(trace)
                trace_net, trace_im, trace_fm = petri.utils.construct_trace_net(trace)
                sync_prod_net, sync_initial_marking, sync_final_marking = petri.synchronous_product.construct(trace_net,
                                                                                                              trace_im,
                                                                                                              trace_fm,
                                                                                                              net,
                                                                                                              i_marking,
                                                                                                              f_marking,
                                                                                                              SKIP)
                reachability_graph = construct_reachability_graph(sync_prod_net, sync_initial_marking)
                print(reachability_graph)
