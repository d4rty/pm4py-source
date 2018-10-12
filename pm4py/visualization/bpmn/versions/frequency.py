from pm4py.visualization.bpmn.util.bpmn_to_figure import bpmn_diagram_to_figure
from pm4py.objects.conversion.bpmn_to_petri import factory as bpmn_converter

def apply(bpmn_graph, parameters=None, bpmn_aggreg_statistics=None):
    """
    Visualize a BPMN graph from a BPMN graph, decorated with frequency, using the given parameters

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
    if parameters is None:
        parameters = {}

    format = parameters["format"] if "format" in parameters else "png"

    file_name = bpmn_diagram_to_figure(bpmn_graph, format, bpmn_aggreg_statistics=bpmn_aggreg_statistics)
    return file_name

def apply_petri(net, initial_marking, final_marking, log=None, bpmn_aggreg_statistics=None, parameters=None):
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
    bpmn_aggreg_statistics
        (Optional) element-wise statistics that should be represented on the BPMN graph
    parameters
        Possible parameters of the algorithm, including:
            format -> Format of the image to render (pdf, png, svg)

    Returns
    -----------
    file_name
        Path of the figure in which the rendered BPMN has been saved
    """

    if parameters is None:
        parameters = {}

    format = parameters["format"] if "format" in parameters else "png"

    bpmn_graph, elements_correspondence, inv_elements_correspondence, el_corr_keys_map = bpmn_converter.apply(net,
                                                                                                              initial_marking,
                                                                                                              final_marking)

    file_name = bpmn_diagram_to_figure(bpmn_graph, format, bpmn_aggreg_statistics=None)
    return file_name