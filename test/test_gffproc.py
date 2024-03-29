import os
import types

import pytest

from tb2neo.gffproc import parse_gff, get_locus_tags, map_functional_category
from test_cli import TEST_GFF, MYCO_TEST_GFF


@pytest.mark.skip(reason="heavy on mem")
def test_parse_gff():
    result = parse_gff(TEST_GFF)
    assert os.path.exists(result) == os.path.exists(TEST_GFF)


def test_get_locus_tags():
    result = get_locus_tags(TEST_GFF, 300)
    assert isinstance(result, types.GeneratorType) is True
    assert len(list(result)) > 0


def test_map_functional_category():
    map_functional_category(gff=MYCO_TEST_GFF)
