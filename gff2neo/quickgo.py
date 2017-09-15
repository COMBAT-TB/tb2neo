"""
Interface to the quickGO interface.
"""
from __future__ import print_function

from bioservices import QuickGO


def fetch_quick_go_data(go_id):
    """
    Retrieve information given a GO identifier.
    :param go_id:
    :return:
    """
    s = QuickGO()
    go_is_a = []
    if not go_id.startswith('GO:'):
        go_id = go_id[go_id.find('[') + 1:]
    result = s.Term(go_id, frmt="obo").split('\n')
    for res in result:
        if 'is_a' in res:
            go_is_a.append(res)
    return go_is_a
