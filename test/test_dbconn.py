from gff2neo.dbconn import graph
from model.core import Gene


def test_db_nodes():
    result = graph.node_labels
    assert "Gene" in result
    assert "Transcript" in result
    assert "CDS" in result


def test_rv0001():
    gene = Gene.select(graph, 'gene:Rv0001').first()
    assert gene.name == 'dnaA'


def test_gene_part_of_transcript():
    part_of = graph.evaluate("MATCH (g:Gene) OPTIONAL MATCH "
                             "((g)<-[:PART_OF]-(t)) RETURN t.parent")
    assert part_of is not None
    located_at = graph.evaluate("MATCH (g:Gene) OPTIONAL MATCH "
                                "((g)-[:LOCATED_AT]->(l)) RETURN l.fmax")
    assert located_at is not None
