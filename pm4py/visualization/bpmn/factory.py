from pm4py.visualization.bpmn.versions import wo_decoration, frequency, performance
from pm4py.visualization.bpmn.util.save_view import *

WO_DECORATION = "wo_decoration"
FREQUENCY_DECORATION = "frequency"
PERFORMANCE_DECORATION = "performance"
FREQUENCY_GREEDY = "frequency_greedy"
PERFORMANCE_GREEDY = "performance_greedy"

VERSIONS = {WO_DECORATION: wo_decoration.apply, FREQUENCY_DECORATION: frequency.apply,
            PERFORMANCE_DECORATION: performance.apply, FREQUENCY_GREEDY: frequency.apply,
            PERFORMANCE_GREEDY: performance.apply}

VERSIONS_PETRI = {WO_DECORATION: wo_decoration.apply_petri, FREQUENCY_DECORATION: frequency.apply_petri, PERFORMANCE_DECORATION: performance.apply_petri,
                  FREQUENCY_GREEDY: frequency.apply_petri_greedy, PERFORMANCE_GREEDY: performance.apply_petri_greedy}

def apply(bpmn_graph, bpmn_aggreg_statistics=None, parameters=None, variant="wo_decoration"):
    """
    Factory method to visualize a BPMN graph using the given parameters

    Parameters
    -----------
    bpmn_graph
        BPMN graph object
    bpmn_aggreg_statistics
        Element-wise statistics that should be represented on the BPMN graph
    parameters
        Possible parameters, of the algorithm, including:
            format -> Format of the image to render (pdf, png, svg)

    Returns
    ----------
    file_name
        Path of the figure in which the rendered BPMN has been saved
    """
    return VERSIONS[variant](bpmn_graph, bpmn_aggreg_statistics=bpmn_aggreg_statistics, parameters=parameters)

def apply_petri(net, initial_marking, final_marking, log=None, aggregated_statistics=None, parameters=None, variant="wo_decoration"):
    """
    Visualize a BPMN graph from a Petri net using the given parameters

    Parameters
    -----------
    net
        Petri net
    initial_marking
        Initial marking
    final_marking
        Final marking
    log
        (Optional) log where the replay technique should be applied
    aggregated_statistics
        (Optional) element-wise statistics calculated on the Petri net
    parameters
        Possible parameters of the algorithm, including:
            format -> Format of the image to render (pdf, png, svg)
            aggregationMeasure -> Measure to use to aggregate statistics
            pm4py.util.constants.PARAMETER_CONSTANT_ACTIVITY_KEY -> Specification of the activity key (if not concept:name)
            pm4py.util.constants.PARAMETER_CONSTANT_TIMESTAMP_KEY -> Specification of the timestamp key (if not time:timestamp)

    Returns
    -----------
    file_name
        Path of the figure in which the rendered BPMN has been saved
    """
    return VERSIONS_PETRI[variant](net, initial_marking, final_marking, log=log, aggregated_statistics=aggregated_statistics, parameters=parameters)