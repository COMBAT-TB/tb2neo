"""
Interface to the quickGO interface.
"""
from __future__ import print_function


def fetch_quick_go_data(quick_go, go_id):
    """
    Retrieve information given a GO identifier.
    :param quick_go:
    :param go_id:
    :return:
    """
    go_is_a = []
    if go_id is None or go_id.startswith("GO:") is False:
        raise ValueError("GO id can't be: {}".format(go_id))
    result = quick_go.Term(go_id, frmt="obo")
    if result:
        for res in result.split('\n'):
            if 'is_a' in res:
                go_is_a.append(res)
    return go_is_a
