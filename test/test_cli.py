import os

from click.testing import CliRunner

from gff2neo.cli import default_gff, check_csv, examine_gff, load_gff


def test_default_gff():
    assert os.stat(default_gff()).st_size > 0


def test_check_csv():
    uniprot_data_csv = "data/uniprot_data.csv"
    assert check_csv(uniprot_data_csv) is True


def test_examine_gff():
    runner = CliRunner()
    result = runner.invoke(examine_gff)
    assert result.exit_code == 0
    assert 'gff_type' in result.output


def test_load_gff():
    runner = CliRunner()
    result = runner.invoke(load_gff)
    assert result.exit_code == 0
