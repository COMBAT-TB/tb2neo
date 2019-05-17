import types

import pytest

from tb2neo.gffproc import get_locus_tags
from tb2neo.uniprot import eu_mapping, query_uniprot
from test_cli import TEST_GFF


def test_search_uniprot():
    identifiers = get_locus_tags(TEST_GFF, 300)
    assert isinstance(identifiers, types.GeneratorType) is True
    result = query_uniprot(locus_tags=identifiers,
                           taxonomy='83332', proteome='UP000001584')
    assert isinstance(result, list) is True


@pytest.mark.parametrize("test_input,expected", [
    (eu_mapping("P13029", to="PDB_ID")[0], "1U2J"),
    (eu_mapping("P9WNW3", to="PDB_ID"), None),
    (eu_mapping("P9WH13", to="TUBERCULIST_ID")[0], "Rv0334"),
    (eu_mapping("P9WPD1", to="TUBERCULIST_ID")[0], "Rv0384c"),

])
def test_eu_mapping(test_input, expected):
    assert test_input == expected
