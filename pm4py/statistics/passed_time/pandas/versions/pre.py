from pm4py.util.constants import PARAMETER_CONSTANT_CASEID_KEY, PARAMETER_CONSTANT_ACTIVITY_KEY, \
    PARAMETER_CONSTANT_TIMESTAMP_KEY
from pm4py.objects.log.util.xes import DEFAULT_NAME_KEY, DEFAULT_TIMESTAMP_KEY
from pm4py.algo.filtering.common.filtering_constants import CASE_CONCEPT_NAME
from pm4py.algo.discovery.dfg.adapters.pandas import df_statistics


def apply(df, activity, parameters=None):
    """
    Gets the time passed from each preceding activity

    Parameters
    -------------
    df
        Dataframe
    activity
        Activity that we are considering
    parameters
        Possible parameters of the algorithm

    Returns
    -------------
    dictio
        Dictionary containing a 'pre' key with the
        list of aggregates times from each preceding activity to the given activity
    """
    if parameters is None:
        parameters = {}

    case_id_glue = parameters[
        PARAMETER_CONSTANT_CASEID_KEY] if PARAMETER_CONSTANT_CASEID_KEY in parameters else CASE_CONCEPT_NAME
    activity_key = parameters[
        PARAMETER_CONSTANT_ACTIVITY_KEY] if PARAMETER_CONSTANT_ACTIVITY_KEY in parameters else DEFAULT_NAME_KEY
    timestamp_key = parameters[
        PARAMETER_CONSTANT_TIMESTAMP_KEY] if PARAMETER_CONSTANT_TIMESTAMP_KEY in parameters else DEFAULT_TIMESTAMP_KEY

    [dfg_frequency, dfg_performance] = df_statistics.get_dfg_graph(df, measure="both", activity_key=activity_key,
                                           case_id_glue=case_id_glue, timestamp_key=timestamp_key)

    pre = []
    sum_perf_pre = 0.0
    sum_acti_pre = 0.0

    for entry in dfg_performance.keys():
        if entry[1] == activity:
            pre.append([entry[0], float(dfg_performance[entry]), int(dfg_frequency[entry])])
            sum_perf_pre = sum_perf_pre + float(dfg_performance[entry]) * float(dfg_frequency[entry])
            sum_acti_pre = sum_acti_pre + float(dfg_frequency[entry])

    perf_acti_pre = 0.0
    if sum_acti_pre > 0:
        perf_acti_pre = sum_perf_pre / sum_acti_pre

    return {"pre": pre, "pre_avg_perf": perf_acti_pre}
