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
    if go_id and go_id.startswith("GO:"):
        result = quick_go.Term(go_id, frmt="obo")
        if result:
            for res in result.split('\n'):
                if 'is_a' in res:
                    go_is_a.append(res)
    else:
        raise ValueError("GO id can't be: {}".format(go_id))
    return go_is_a
