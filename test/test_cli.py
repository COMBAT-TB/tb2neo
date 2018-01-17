import os

import pytest
from click.testing import CliRunner

from gff2neo.cli import default_gff, check_csv, examine_gff, load_gff, \
    load_uniprot_data


@pytest.fixture(scope="module")
def cli_runner():
    runner = CliRunner()
    return runner


def test_default_gff():
    """
    Test if empty
    :return:
    """
    assert os.stat(default_gff()).st_size > 0


def test_check_csv():
    uniprot_data_csv = "data/uniprot_data.csv"
    assert check_csv(uniprot_data_csv) is True


def test_examine_gff(cli_runner):
    result = cli_runner.invoke(examine_gff)
    assert result.exit_code == 0
    assert 'gff_type' in result.output


def test_load_gff(cli_runner):
    result = cli_runner.invoke(load_gff)
    assert result.exit_code == 0


def test_load_uniprot_data(cli_runner):
    result = cli_runner.invoke(load_uniprot_data)
    assert result.exit_code == 0
