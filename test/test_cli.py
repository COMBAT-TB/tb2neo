import os

import pytest
from click.testing import CliRunner

from gff2neo.cli import examine_gff, load_gff, sample_gff, load_uniprot_data, load_publications, load_reactome_pathways, \
    load_go_terms, load_drugbank_data, load_kegg_pathways

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


def test_examine_gff(cli_runner):
    result = cli_runner.invoke(examine_gff, ["--gff_files", sample_gff_dir])
    assert result.exit_code == 0


def test_load_gff(cli_runner):
    result = cli_runner.invoke(load_gff, ["--gff_files", sample_gff_dir])
    assert result.exit_code == 0


def test_load_uniprot_data(cli_runner):
    result = cli_runner.invoke(load_uniprot_data, ["--gff_files", sample_gff_dir])
    assert result.exit_code == 0


def test_load_go_terms(cli_runner):
    result = cli_runner.invoke(load_go_terms)
    assert result.exit_code == 0


def test_load_drugbank_data(cli_runner):
    result = cli_runner.invoke(load_drugbank_data)
    assert result.exit_code == 0


def test_load_kegg_pathways(cli_runner):
    result = cli_runner.invoke(load_kegg_pathways)
    assert result.exit_code == 0


def test_load_reactome_pathways(cli_runner):
    result = cli_runner.invoke(load_reactome_pathways)
    assert result.exit_code == 0


def test_load_publications(cli_runner):
    result = cli_runner.invoke(load_publications)
    assert result.exit_code == 0
