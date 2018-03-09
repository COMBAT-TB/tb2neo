import pytest

from gff2neo.dbconn import graph, split_gene_names
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
    ("MATCH (g:Gene) OPTIONAL MATCH ((g)<-[:PART_OF]-(t)) RETURN t.parent", type(None)),
    ("MATCH (g:Gene) OPTIONAL MATCH ((g)-[:LOCATED_AT]->(l)) RETURN l.fmax", type(None)),
    ("MATCH (g:Gene) OPTIONAL MATCH ((g)-[:ENCODES]->(p)) RETURN p.parent", type(None)),
])
def test_db_data(test_input, expected):
    assert isinstance(type(graph.evaluate(test_input)), expected) is False


@pytest.mark.parametrize("test_input,expected", [
    (split_gene_names("MT0511/MT0512"), list),
    (split_gene_names("MT1076 MT1237 MT3197"), list),
    (split_gene_names("MT0511;MT0512"), list),
    (split_gene_names(None), type(None)),
])
def test_split_gene_name(test_input, expected):
    assert isinstance(test_input, expected) is True
