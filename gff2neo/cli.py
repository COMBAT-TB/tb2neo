import click

from gffproc import *


def default_gff():
    """
    Check if gff is given, use default if not.
    :return:
    """
    gff_file = "data/sample.gff3"
    if os.path.isdir(os.getcwd() + "/data") and os.path.exists(gff_file):
        click.echo("\nUsing {}.\n".format(gff_file), nl=True, err=True)
    else:
        raise OSError("Couldn't find GFF file: {}!".format(gff_file))
    return gff_file


def check_csv(csvfile):
    """
    Check if csv file exists and is not empty
    :param csvfile:
    :return:
    """
    _file = False
    if os.path.exists(csvfile) and os.stat(csvfile).st_size > 0:
        _file = True
    return _file


@click.group()
def cli():
    """
    This script parses a GFF file and builds a Neo4j Graph database.
    """
    pass


@cli.command()
@click.argument('gff_file', required=False, type=click.Path(exists=True, file_okay=True))
def examine_gff(gff_file):
    """
    Examine GFF file features.
    :param gff_file
    :return
    """
    if gff_file is None:
        gff_file = default_gff()
    examine_gff_file(gff_file=gff_file)


@cli.command()
@click.argument('gff_file', required=False, type=click.Path(exists=True, file_okay=True))
def load_gff(gff_file):
    """
    Load GFF features and build relationships.
    :param gff_file:
    :return:
    """
    # Deleting existing data
    delete_db_data()
    if gff_file is None:
        gff_file = default_gff()
    parse_gff(gff_file)
    build_gff_relationships()


@cli.command()
@click.argument('gff_file', required=False, type=click.Path(exists=True, file_okay=True))
def load_uniprot_data(gff_file):
    """
    Load UniProt data using GFF RVs (locus tags).
    :param gff_file:
    :return:
    """
    if check_csv(uniprot_data_csv):
        create_protein_nodes()
    else:
        if gff_file is None:
            gff_file = default_gff()
        query_uniprot(get_locus_tags(gff_file=gff_file, chunk=400))
        create_protein_nodes()
        map_gene_to_protein(get_locus_tags(gff_file, 400))


@cli.command()
@click.argument('gff_file', required=False, type=click.Path(exists=True, file_okay=True))
def load_drugbank_data(gff_file):
    """
    Load DrugBank data.
    :param gff_file:
    :return:
    """
    # Let's check if we have Proteins to map to.
    proteins = graph.data("MATCH (p:Protein) RETURN p.entry_name LIMIT 1")
    if len(proteins) > 0 and check_csv(uniprot_data_csv) \
            and check_csv(target_protein_ids_csv) and check_csv(drug_vocab_csv):
        create_drugbank_nodes()
    else:
        if gff_file is None:
            gff_file = default_gff()
        query_uniprot(get_locus_tags(gff_file=gff_file, chunk=400))
        create_protein_nodes()
        map_gene_to_protein(get_locus_tags(gff_file, 400))
        create_drugbank_nodes()


@cli.command()
@click.argument('gff_file', required=False, type=click.Path(exists=True, file_okay=True))
def load_go_terms(gff_file):
    """
    Load GO terms.
    :param gff_file:
    :return:
    """
    if check_csv(uniprot_data_csv):
        create_go_term_nodes()
    else:
        if gff_file is None:
            gff_file = default_gff()
        query_uniprot(get_locus_tags(gff_file=gff_file, chunk=400))
        create_go_term_nodes()


@cli.command()
@click.argument('gff_file', required=False, type=click.Path(exists=True, file_okay=True))
def load_publications(gff_file):
    """
    Load Publications.
    :param gff_file:
    :return:
    """
    if check_csv(uniprot_data_csv):
        create_publication_nodes()
    else:
        if gff_file is None:
            gff_file = default_gff()
        query_uniprot(get_locus_tags(gff_file=gff_file, chunk=400))
        create_publication_nodes()


@cli.command()
@click.argument('gff_file', required=False, type=click.Path(exists=True, file_okay=True))
def load_kegg_pathways(gff_file):
    """
    Load KEGG Pathways
    :param gff_file:
    :return:
    """
    # Let's check if we have Proteins to map to.
    proteins = graph.data("MATCH (p:Protein) RETURN p.entry_name LIMIT 1")
    if len(proteins) > 0:
        create_kegg_pathways_nodes()
    else:
        if gff_file is None:
            gff_file = default_gff()
        query_uniprot(get_locus_tags(gff_file=gff_file, chunk=400))
        create_kegg_pathways_nodes()


@cli.command()
@click.argument('gff_file', required=False, type=click.Path(exists=True, file_okay=True))
def load_reactome_pathways(gff_file):
    """
    Load REACTOME Pathways
    :param gff_file:
    :return:
    """
    # Let's check if we have Proteins to map to.
    proteins = graph.data("MATCH (p:Protein) RETURN p.entry_name LIMIT 1")
    if len(proteins) > 0 and check_csv(uniprot_data_csv):
        create_reactome_pathway_nodes()
    else:
        if gff_file is None:
            gff_file = default_gff()
        query_uniprot(get_locus_tags(gff_file=gff_file, chunk=400))
        create_reactome_pathway_nodes()


if __name__ == '__main__':
    cli()
