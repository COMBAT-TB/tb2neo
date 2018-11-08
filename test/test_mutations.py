import pytest

from gff2neo.mutations import get_drug_info


@pytest.mark.parametrize("test_input,expected", [
    (get_drug_info("amikacin"), "DB00479"),
    (get_drug_info("ethionamide"), "DB00609"),
    (get_drug_info("amikacin"), "DB00479"),
])
def test_get_drug_info(test_input, expected):
    assert test_input == expected
