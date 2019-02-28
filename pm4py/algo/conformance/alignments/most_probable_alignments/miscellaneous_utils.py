import math

from pm4py.algo.conformance.alignments.utils import SKIP


def apply_log_transformation(probability):
    if probability is None:
        raise ValueError('Probability is of type NoneType')
    # TODO optimize --> apply log transformation not on the fly but apply it if possible before
    # --> try to avoid duplicate calls
    if not (0 <= probability <= 1):
        raise ValueError('probability needs to be in the range [0,1]')
    if probability > 0:
        return - math.log(probability)
    elif probability == 0:
        # -log(0) ~ math.inf
        return math.inf
    else:
        raise ValueError('Log Transformation not applicable to numbers below 0')


def sync_product_net_place_belongs_to_process_net(p):
    """
    :param p: Place object; represents a place in a synchronous product net
    :return: Boolean - true if the synchronous product net state represents a state of the process net
    """
    return p.name[0] == SKIP and p.name[1] != SKIP


def is_model_move(t, skip):
    return t.label[0] == skip and t.label[1] != skip


def is_log_move(t, skip):
    return t.label[0] != skip and t.label[1] == skip