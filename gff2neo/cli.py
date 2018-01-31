import click

from gff2neo.gffproc import *

sample_gff = os.getcwd() + "/data/sample_gff/sample.gff3"


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
    'sample': {'taxonomy': '83332', 'proteome': 'UP000001584'},
}


@click.group()
def cli():
    """
    This script parses a GFF file and builds a Neo4j Graph database.
    """
    pass


@cli.command()
@click.option('gff_files', '--gff_files', type=click.Path(exists=True), default=sample_gff,
              help="Path to GFF file(s) (default: '{}')".format(sample_gff))
def examine_gff(gff_files):
    """
    Examine features from GFF file.
    :param gff_files: GFF file directory
    :return
    """
    if not os.path.isdir(gff_files):
        gff_files = os.path.dirname(sample_gff)
    for gff_file in os.listdir(os.path.abspath(gff_files)):
        if gff_file.endswith(".gff3"):
            click.secho("Examining: {}".format(gff_file), fg="green")
            examine_gff_file(gff_file=os.path.abspath(gff_files) + "/" + gff_file)


@cli.command()
@click.option('gff_files', '--gff_files', type=click.Path(exists=True), default=sample_gff,
              help="Path to GFF file(s) (default: '{}')".format(sample_gff))
def load_gff(gff_files):
    """
    Load GFF features and build relationships.
    :param gff_files:
    :return:
    """
    # Deleting existing data
    delete_db_data()
    if not os.path.isdir(gff_files):
        gff_files = os.path.dirname(sample_gff)
    for gff_file in os.listdir(os.path.abspath(gff_files)):
        if gff_file.endswith(".gff3"):
            click.secho("Loading {}".format(os.path.abspath(gff_files) + "/" + gff_file), fg="green")
            parse_gff(os.path.abspath(gff_files) + "/" + gff_file)
            build_gff_relationships()


@cli.command()
@click.option('gff_files', '--gff_files', type=click.Path(exists=True), default=sample_gff,
              help="Path to GFF file(s) (default: '{}')".format(sample_gff))
def load_uniprot_data(gff_files):
    """
    Load UniProt data using GFF (locus tags).
    :param gff_files:
    :return:
    """
    taxonomy = None
    proteome = None
    click.secho("Loading UniProt data...", fg="green")
    if check_csv(uniprot_data_csv):
        click.secho("Found CSV data...", fg="green")
        create_protein_nodes()
    else:
        if not os.path.isdir(gff_files):
            gff_files = os.path.dirname(sample_gff)
        for gff_file in os.listdir(os.path.abspath(gff_files)):
            gff_file = os.path.abspath(gff_files) + "/" + gff_file
            if gff_file.endswith(".gff3"):
                for k, v in uniprot_config.items():
                    if k in gff_file.split("/")[-1]:
                        taxonomy = v['taxonomy']
                        proteome = v['proteome']
                click.secho("Fetchig data about {}".format(gff_file), fg="green")
                locus_tags = get_locus_tags(gff_file=gff_file, chunk=400)
                query_uniprot(locus_tags=locus_tags, taxonomy=taxonomy, proteome=proteome)
                # TODO: Need to refactor
                create_protein_nodes()
                map_gene_to_protein(get_locus_tags(gff_file, 400))


@cli.command()
@click.option('gff_files', '--gff_files', type=click.Path(exists=True), default=sample_gff,
              help="Path to GFF file(s) (default: '{}')".format(sample_gff))
def load_drugbank_data(gff_files):
    """
    Load DrugBank data.
    :param gff_files:
    :return:
    """
    # Let's check if we have Proteins to map to.
    proteins = graph.data("MATCH (p:Protein) RETURN p.entry_name LIMIT 1")
    if len(proteins) > 0 and check_csv(uniprot_data_csv) \
            and check_csv(target_protein_ids_csv) and check_csv(drug_vocab_csv):
        create_drugbank_nodes()
    else:
        if not os.path.isdir(gff_files):
            gff_files = os.path.dirname(sample_gff)
        for gff_file in gff_files:
            query_uniprot(get_locus_tags(gff_file=gff_file, chunk=400))
        create_protein_nodes()
        for gff_file in gff_files:
            map_gene_to_protein(get_locus_tags(gff_file, 400))
        create_drugbank_nodes()


@cli.command()
@click.option('gff_files', '--gff_files', type=click.Path(exists=True), default=sample_gff,
              help="Path to GFF file(s) (default: '{}')".format(sample_gff))
def load_go_terms(gff_files):
    """
    Load GO terms.
    :param gff_files:
    :return:
    """
    if check_csv(uniprot_data_csv):
        create_go_term_nodes()
    else:
        if not os.path.isdir(gff_files):
            gff_files = os.path.dirname(sample_gff)
        for gff_file in gff_files:
            query_uniprot(get_locus_tags(gff_file=gff_file, chunk=400))
        create_go_term_nodes()


@cli.command()
@click.option('gff_files', '--gff_files', type=click.Path(exists=True), default=sample_gff,
              help="Path to GFF file(s) (default: '{}')".format(sample_gff))
def load_publications(gff_files):
    """
    Load Publications.
    :param gff_files:
    :return:
    """
    if check_csv(uniprot_data_csv):
        create_publication_nodes()
    else:
        if not os.path.isdir(gff_files):
            gff_files = os.path.dirname(sample_gff)
        for gff_file in gff_files:
            query_uniprot(get_locus_tags(gff_file=gff_file, chunk=400))
        create_publication_nodes()


@cli.command()
@click.option('gff_files', '--gff_files', type=click.Path(exists=True), default=sample_gff,
              help="Path to GFF file(s) (default: '{}')".format(sample_gff))
def load_kegg_pathways(gff_files):
    """
    Load KEGG Pathways
    :param gff_files:
    :return:
    """
    # Let's check if we have Proteins to map to.
    proteins = graph.data("MATCH (p:Protein) RETURN p.entry_name LIMIT 1")
    if len(proteins) > 0:
        create_kegg_pathways_nodes()
    else:
        if not os.path.isdir(gff_files):
            gff_files = os.path.dirname(sample_gff)
        for gff_file in gff_files:
            query_uniprot(get_locus_tags(gff_file=gff_file, chunk=400))
        create_kegg_pathways_nodes()


@cli.command()
@click.option('gff_files', '--gff_files', type=click.Path(exists=True), default=sample_gff,
              help="Path to GFF file(s) (default: '{}')".format(sample_gff))
def load_reactome_pathways(gff_files):
    """
    Load REACTOME Pathways
    :param gff_files:
    :return:
    """
    # Let's check if we have Proteins to map to.
    proteins = graph.data("MATCH (p:Protein) RETURN p.entry_name LIMIT 1")
    if len(proteins) > 0 and check_csv(uniprot_data_csv):
        create_reactome_pathway_nodes()
    else:
        if not os.path.isdir(gff_files):
            gff_files = os.path.dirname(sample_gff)
        for gff_file in gff_files:
            query_uniprot(get_locus_tags(gff_file=gff_file, chunk=400))
        create_reactome_pathway_nodes()


if __name__ == '__main__':
    cli()
