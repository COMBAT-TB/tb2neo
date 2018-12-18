"""
Testing ncbi module
"""
import types

import pytest

from gff2neo.ncbi import fetch_publication_list, get_fasta, search_pubmed

pubmed_ids = "8733228; 9634230; 21969609; 21219854; 22545130; 26045430"
pubmed_id_list = [p for p in pubmed_ids.split("; ") if p is not '']


def test_search_pubmed():
    res = search_pubmed(genename="rmlA")
    assert isinstance(res, list)
    assert len(res) is 10


@pytest.mark.parametrize("test_input,expected", [
    (fetch_publication_list(pubmed_id_list), types.GeneratorType),
    (fetch_publication_list(pubmed_id_list, rettype='xml'),
     types.GeneratorType),
])
def test_fetch_publication_list(test_input, expected):
    assert isinstance(test_input, expected) is True


def test_get_fasta():
    res = get_fasta(strain="h37rv")
    assert res[:4] == 'TTGA'
