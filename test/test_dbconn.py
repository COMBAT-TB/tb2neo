"""
Test dbconn module
"""
import pytest

from tb2neo.dbconn import (create_chromosome_nodes, create_publication_nodes,
                           graph, split_gene_names)
from combattbmodel.core import Chromosome, Gene
from test_cli import UNIPROT_DATA


@pytest.mark.skip(reason="heavy on mem")
def test_create_chromosome_nodes():
    create_chromosome_nodes(strain="h37rv")
    chromosome = Chromosome.select(graph).first()
    assert chromosome.residues[:4] == "TTGA"


def test_db_nodes():
    result = graph.node_labels
    assert "Organism" in result
    assert "Chromosome" in result
    assert "Gene" in result
    assert "MRna" in result
    assert "CDS" in result
    assert "Protein" in result


def test_rv0001():
    gene = Gene.select(graph, "Rv0001").first()
    assert gene.name == "dnaA"
    assert gene.category != ""
    assert gene.residues != ""


@pytest.mark.parametrize("test_input,expected", [
    ("MATCH (g:Gene) OPTIONAL MATCH ((g)<-[:PART_OF]-(t)) RETURN t.parent",
     type(None)),
    ("MATCH (g:Gene) OPTIONAL MATCH ((g)-[:LOCATED_AT]->(l)) RETURN l.fmax",
     type(None)),
    ("MATCH (g:Gene) OPTIONAL MATCH ((g)-[:ENCODES]->(p)) RETURN p.parent",
     type(None)),
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


def test_create_publication_nodes():
    res = create_publication_nodes(uniprot_data=UNIPROT_DATA)
    assert isinstance(res, set)


# @pytest.mark.skip(reason="")
def test_ppi_score():
    cypher_q = f"MATCH (pa:Protein {{ uniquename: 'P9WQD9' }})"
    cypher_q += f"-[r:INTERACTS_WITH]->"
    cypher_q += f"(pb:Protein {{ uniquename: 'P9WIE5' }}) "
    cypher_q += f"RETURN r.score"
    data = graph.run(cypher_q).data()[0]
    assert data['r.score'] == 0.881
