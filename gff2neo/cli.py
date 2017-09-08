import click

from gffproc import *


def default_gff():
    gff_file_ = "data/MTB_H37rv.gff3"
    if os.path.isdir(os.getcwd() + "/data") and os.path.exists(gff_file_):
        click.echo("\nUsing {}.\n".format(gff_file_), nl=True, err=True)
    else:
        raise OSError("Couldn't find GFF file: {}!".format(gff_file_))
    return gff_file_


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
    examine(gff_file=gff_file)


@cli.command()
@click.argument('gff_file', required=False, type=click.Path(exists=True, file_okay=True))
def load_gff(gff_file):
    """
    Load GFF features and build relationships.
    :param gff_file:
    :return:
    """
    # Deleting existing data
    delete_data()
    if gff_file is None:
        gff_file = default_gff()
    parse_gff(gff_file)
    build_gff_rels()


@cli.command()
@click.argument('gff_file', required=False, type=click.Path(exists=True, file_okay=True))
def load_uniprot_data(gff_file):
    """
    Load UniProt data using GFF RVs (locus tags).
    :param gff_file:
    :return:
    """
    if os.path.exists(uniprot_data_csv) and os.stat(uniprot_data_csv).st_size > 0:
        create_uniprot_nodes()
    else:
        if gff_file is None:
            gff_file = default_gff()
        query_uniprot(get_locus_tags(gff_file=gff_file, chunk=400))
        create_uniprot_nodes()


@cli.command()
@click.argument('gff_file', required=False, type=click.Path(exists=True, file_okay=True))
def load_go_terms(gff_file):
    """
    Load GO terms.
    :param gff_file:
    :return:
    """
    if os.path.exists(uniprot_data_csv) and os.stat(uniprot_data_csv).st_size > 0:
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
    if os.path.exists(uniprot_data_csv) and os.stat(uniprot_data_csv).st_size > 0:
        create_pub_nodes()
    else:
        if gff_file is None:
            gff_file = default_gff()
        query_uniprot(get_locus_tags(gff_file=gff_file, chunk=400))
        create_pub_nodes()


if __name__ == '__main__':
    cli()
