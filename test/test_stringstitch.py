import pytest

from tb2neo.stringstitch import fetch_string_data


@pytest.mark.parametrize("test_input,expected", [
    (fetch_string_data('Rv0678')[0].get('score'), 0.837),
    (fetch_string_data('P0CW33'), None),
    (fetch_string_data(None)[0].get('score'), 0.481),
])
def test_fetch_string_data(test_input, expected):
    assert test_input == expected
