import types

import pytest

from gff2neo.gffproc import get_locus_tags
from gff2neo.uniprot import query_uniprot, eu_mapping
from .test_cli import TEST_GFF


def test_search_uniprot():
    identifiers = get_locus_tags(TEST_GFF, 300)
    assert isinstance(identifiers, types.GeneratorType) is True
    result = query_uniprot(locus_tags=identifiers, taxonomy='83332', proteome='UP000001584')
    assert isinstance(result, list) is True


@pytest.mark.parametrize("test_input,expected", [
    (eu_mapping("P13029", to="PDB_ID"), list),
    (eu_mapping("P9WNW3", to="PDB_ID"), type(None)),
])
def test_eu_mapping(test_input, expected):
    assert isinstance(test_input, expected) is True


def test_eu_mapping_error():
    with pytest.raises(ValueError):
        eu_mapping(None, to="PDB_ID")
