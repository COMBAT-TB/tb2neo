import click

from gff2neo.gffproc import *
from gff2neo.uniprot import UNIPROT_DATA

CURR_DIR = os.path.dirname(os.path.abspath(__file__))


# UNIPROT_DATA = os.path.join(CURR_DIR, "data/uniprot/uniprot_data.csv")


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


uniprot_config = {
    'h37rv': {'taxonomy': '83332', 'proteome': 'UP000001584'},
    'cdc1551': {'taxonomy': '83331', 'proteome': 'UP000001020'},
}


def get_taxonomy_and_proteome(gff_file):
    """
    Get Taxonomy and Proteome
    :param gff_file:
    :return:
    """
    result = dict()
    file_name = gff_file.split("/")[-1]
    for k, v in uniprot_config.items():
        if k in file_name:
            result['taxonomy'] = v['taxonomy']
            result['proteome'] = v['proteome']
    return result


@click.group()
def cli():
    """
    This script parses a GFF file and builds a Neo4j Graph database.
    """
    pass


@cli.command()
@click.argument('gff_files', type=click.Path(exists=True))
def examine_gff(gff_files):
    """
    Examine features from GFF file.
    :param gff_files: GFF file directory
    :return
    """
    if os.path.isdir(gff_files):
        for root, dirs, files in os.walk(gff_files):
            for gff_file in files:
                gff_file = '/'.join([os.path.abspath(gff_files), gff_file])
                if gff_file.endswith(".gff3"):
                    click.secho("Examining: {}".format(gff_file), fg="green")
                    examine_gff_file(gff_file=gff_file)


@cli.command()
@click.argument('gff_files', type=click.Path(exists=True))
def load_gff(gff_files):
    """
    Load GFF features and build relationships.
    :param gff_files:
    :return:
    """
    # Deleting existing data
    delete_db_data()
    if os.path.isdir(gff_files):
        for root, dirs, files in os.walk(gff_files):
            for gff_file in files:
                gff_file = '/'.join([os.path.abspath(gff_files), gff_file])
                if gff_file.endswith(".gff3"):
                    click.secho("Loading {}".format(gff_file), fg="green")
                    parse_gff(gff_file)
                    build_gff_relationships()


@cli.command()
@click.argument('gff_files', type=click.Path(exists=True))
def load_uniprot_data(gff_files):
    """
    Load UniProt data using GFF (locus tags).
    :param gff_files:
    :return:
    """
    click.secho("\nLoading UniProt data...", fg="green")
    if check_csv(UNIPROT_DATA) and not gff_files:
        click.secho("Found CSV data...", fg="green")
        create_protein_nodes()
    else:
        if os.path.isdir(gff_files):
            for root, dirs, files in os.walk(gff_files):
                for gff_file in files:
                    gff_file = '/'.join([os.path.abspath(gff_files), gff_file])
                    if gff_file.endswith(".gff3"):
                        result = get_taxonomy_and_proteome(gff_file)
                        locus_tags = get_locus_tags(gff_file=gff_file, chunk=400)
                        map_gene_to_orthologs(get_locus_tags(gff_file, 400))
                        query_uniprot(locus_tags=locus_tags, taxonomy=result['taxonomy'], proteome=result['proteome'])
                        # TODO: Need to refactor
                        create_protein_nodes()


@cli.command()
def load_drugbank_data():
    """
    Load DrugBank data.
    :return:
    """
    # Let's check if we have Proteins to map to.
    proteins = graph.data("MATCH (p:Protein) RETURN p.entry_name LIMIT 1")
    if len(proteins) > 0 and check_csv(UNIPROT_DATA) \
            and check_csv(TARGET_PROTEIN_IDS) and check_csv(DRUG_VOCAB):
        create_drugbank_nodes()
    else:
        sys.stderr.write(
            "Unable to load Drugbank data!\n Check if CSV files are in place and that we have a database with Proteins.")


@cli.command()
def load_go_terms():
    """
    Load GO terms.
    :return:
    """
    if check_csv(UNIPROT_DATA):
        create_go_term_nodes()
    else:
        sys.stderr.write(
            "Unable to create GOTerms !\n Check if CSV files are in place.")


@cli.command()
def load_publications():
    """
    Load Publications.
    :return:
    """
    if check_csv(UNIPROT_DATA):
        create_publication_nodes()
    else:
        sys.stderr.write(
            "Unable to load Publications!\n Check if CSV files are in place")


@cli.command()
def load_kegg_pathways():
    """
    Load KEGG Pathways
    :return:
    """
    # Let's check if we have Proteins to map to.
    proteins = graph.data("MATCH (p:Protein) RETURN p.entry_name LIMIT 1")
    if len(proteins) > 0:
        create_kegg_pathways_nodes()
    else:
        sys.stderr.write(
            "Unable to load KEGG Pathways!\n Check if we have a database with Proteins.")


@cli.command()
def load_reactome_pathways():
    """
    Load REACTOME Pathways
    :return:
    """
    # Let's check if we have Proteins to map to.
    proteins = graph.data("MATCH (p:Protein) RETURN p.entry_name LIMIT 1")
    if len(proteins) > 0 and check_csv(UNIPROT_DATA):
        create_reactome_pathway_nodes()
    else:
        sys.stderr.write(
            "Unable to load REACTOME Pathways!\n Check if we have a database with Proteins.")


if __name__ == '__main__':
    cli()
