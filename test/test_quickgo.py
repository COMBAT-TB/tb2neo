"""
Testing quickgo module
"""
import pytest

from tb2neo.quickgo import query_quickgo


# def test_fetch_quick_go_id():
#     with pytest.raises(ValueError):
#         fetch_quick_go_data(quick_go, None)
#         fetch_quick_go_data(quick_go, "0005886")
#
#
# def test_fetch_quick_go_data():
#     result = fetch_quick_go_data(quick_go, "GO:0005886")
#     assert isinstance(result, list) is True

@pytest.mark.parametrize("test_input,expected", [
    (query_quickgo("GO:0005576,GO:0005618,GO:0016021").status_code, 200),
    (query_quickgo("GO:0003993,GO:0005886,GO:0046872").status_code, 200),
    (query_quickgo("GO:0003824,GO:0005886,GO:0016021,GO:0046872,GO:0051539").status_code, 200),
    (query_quickgo(
        "GO:0004160,GO:0005886,GO:0009082,GO:0009097,GO:0009099,GO:0040007,GO:0046872,GO:0051539").status_code, 200),

])
def test_query_quickgo(test_input, expected):
    assert test_input == expected
    assert isinstance(expected, int) is True
