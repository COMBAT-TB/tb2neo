import click

from gffproc import *


@click.group()
def cli():
    """
    This script parses a GFF file and builds a Neo4j Graph database.
    """
    pass


@cli.command()
@click.argument('gff_file', type=click.Path(exists=True, file_okay=True))
def inspect_gff(gff_file):
    """Examine GFF File."""
    examine(gff_file=gff_file)


if __name__ == '__main__':
    cli()
