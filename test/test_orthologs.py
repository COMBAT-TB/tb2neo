import pytest

from gff2neo.orthologs import fetch_ortholog


def test_fetch_ortholog():
    with pytest.raises(ValueError):
        fetch_ortholog()
    with pytest.raises(Exception):
        # Looks like the service is down!
        fetch_ortholog(locus_tag="Rv0239")
