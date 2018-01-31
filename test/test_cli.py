import os

import pytest
from click.testing import CliRunner

from gff2neo.cli import examine_gff, load_gff, sample_gff, load_uniprot_data

sample_gff_dir = os.path.dirname(sample_gff)


@pytest.fixture(scope="module")
def cli_runner():
    runner = CliRunner()
    return runner


def test_default_gff():
    """
    Test if empty
    :return:
    """
    assert os.stat(sample_gff).st_size > 0


@pytest.mark.parametrize("test_input,expected_output", [
    (examine_gff, 0),
])
def test_examine_gff(cli_runner, test_input, expected_output):
    result = cli_runner.invoke(test_input)
    assert result.exit_code == expected_output


@pytest.mark.parametrize("test_input,expected_output", [
    (load_gff, 0),
])
def test_load_gff(cli_runner, test_input, expected_output):
    result = cli_runner.invoke(test_input)
    assert result.exit_code == expected_output


def test_load_uniprot_data(cli_runner):
    result = cli_runner.invoke(load_uniprot_data, ["--gff_files", os.getcwd() + "/data/sample_gff/"])
    assert result.exit_code == 0

# def test_load_go_terms(cli_runner):
#     result = cli_runner.invoke(load_go_terms)
#     assert result.exit_code == 0
#
#
# def test_load_drugbank_data(cli_runner):
#     result = cli_runner.invoke(load_drugbank_data)
#     assert result.exit_code == 0
#
#
# def test_load_reactome_pathways(cli_runner):
#     result = cli_runner.invoke(load_reactome_pathways)
#     assert result.exit_code == 0
