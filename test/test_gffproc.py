import os
import types

from gff2neo.cli import sample_gff
from gff2neo.gffproc import parse_gff, get_locus_tags


def test_parse_gff():
    result = parse_gff(sample_gff)
    assert os.path.exists(result) == os.path.exists(sample_gff)


def test_get_locus_tags():
    result = get_locus_tags(sample_gff, 300)
    assert isinstance(result, types.GeneratorType) is True
    assert len(list(result)) > 0
