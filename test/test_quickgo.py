import pytest

from gff2neo.cli import quick_go
from gff2neo.quickgo import fetch_quick_go_data


def test_fetch_quick_go_data():
    result = fetch_quick_go_data(quick_go, "GO:0005886")
    assert isinstance(result, list) is True
    with pytest.raises(ValueError):
        fetch_quick_go_data(quick_go, None)
        fetch_quick_go_data(quick_go, "0005886")
