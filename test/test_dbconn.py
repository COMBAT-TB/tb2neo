import types

import pytest

from gff2neo.dbconn import graph
from model.core import Gene


def test_db_nodes():
    result = graph.node_labels
    assert "Gene" in result
    assert "Transcript" in result
    assert "CDS" in result


def test_rv0001():
    gene = Gene.select(graph, 'Rv0001').first()
    assert gene.name == 'dnaA'


@pytest.mark.parametrize("test_input,expected", [
    ("MATCH (g:Gene) OPTIONAL MATCH ((g)<-[:PART_OF]-(t)) RETURN t.parent", types.NoneType),
    ("MATCH (g:Gene) OPTIONAL MATCH ((g)-[:LOCATED_AT]->(l)) RETURN l.fmax", types.NoneType),
    ("MATCH (g:Gene) OPTIONAL MATCH ((g)-[:ENCODES]->(p)) RETURN p.parent", types.NoneType),
])
def test_db_data(test_input, expected):
    assert isinstance(type(graph.evaluate(test_input)), expected) is False
