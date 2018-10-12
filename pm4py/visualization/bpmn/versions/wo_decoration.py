from pm4py.visualization.bpmn.util.bpmn_to_figure import bpmn_diagram_to_figure

def apply(bpmn_graph, parameters=None):
    """
    Visualize a BPMN graph from a BPMN graph using the given parameters

    Parameters
    -----------
    bpmn_graph
        BPMN graph object
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

    file_name = bpmn_diagram_to_figure(bpmn_graph, format, bpmn_aggreg_statistics=None)
    return file_name