"""
Testing CLI module
"""
import os

import pytest
from click.testing import CliRunner

from gff2neo.cli import delete, examine_gff, load_gff, load_uniprot_data, load_organism, \
    load_publications, load_reactome_pathways, load_go_terms, \
    load_drugbank_data, load_kegg_pathways, load_chromosome, \
    load_known_mutations, load_operons, load_srna_data, SRNA_TXTFILE


CURR_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_GFF = os.path.join(CURR_DIR, "test_gff/h37rv-sample.gff3")
MYCO_TEST_GFF = os.path.join(CURR_DIR, "test_myco_gff/test_myco_h37rv.gff")
TEST_GFF_DIR = os.path.join(CURR_DIR, "test_gff/")
TEST_OPERON_FILE = os.path.join(
    CURR_DIR, "test/test_operon/h37rv_operon_sample.txt")
TEST_OPERON_FILE_DIR = os.path.join(CURR_DIR, "test_operon/")
TEST_MUTATIONS_DIR = os.path.join(CURR_DIR, "test_mutations/")
UNIPROT_DATA = os.path.join(CURR_DIR, "test_uniprot_data/uniprot_data.csv")


@pytest.fixture(scope="module")
def cli_runner():
    runner = CliRunner()
    return runner


def test_default_gff():
    """
    Test if empty
    :return:
    """
    assert os.stat(TEST_GFF).st_size > 0


def test_delete(cli_runner):
    result = cli_runner.invoke(delete)
    assert result.exit_code == 0


def test_load_organism(cli_runner):
    result = cli_runner.invoke(load_organism, [TEST_GFF_DIR])
    assert result.exit_code == 0

# @pytest.mark.skip(reason="no way of currently testing this")


def test_load_chromosome(cli_runner):
    result = cli_runner.invoke(load_chromosome, ["h37rv"])
    assert result.exit_code == 0


def test_examine_gff(cli_runner):
    result = cli_runner.invoke(examine_gff, [TEST_GFF_DIR])
    assert result.exit_code == 0


def test_load_gff(cli_runner):
    result = cli_runner.invoke(load_gff, [TEST_GFF_DIR])
    assert result.exit_code == 0


@pytest.mark.skip(reason="no way of currently testing this")
def test_load_operons(cli_runner):
    result = cli_runner.invoke(load_operons, [TEST_OPERON_FILE_DIR])
    assert result.exit_code == 0


@pytest.mark.skip(reason="no way of currently testing this")
def test_load_known_mutations(cli_runner):
    result = cli_runner.invoke(load_known_mutations, [TEST_MUTATIONS_DIR])
    assert result.exit_code == 0


def test_load_uniprot_data(cli_runner):
    result = cli_runner.invoke(load_uniprot_data, [TEST_GFF_DIR])
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


@pytest.mark.skip(reason="'NoneType' object is not subscriptable")
def test_load_reactome_pathways(cli_runner):
    result = cli_runner.invoke(load_reactome_pathways)
    assert result.exit_code == 0


def test_load_publications(cli_runner):
    result = cli_runner.invoke(load_publications)
    assert result.exit_code == 0


def test_load_srna_data(cli_runner):
    result = cli_runner.invoke(load_srna_data, [SRNA_TXTFILE])
    assert result.exit_code == 0
