from setuptools import setup

setup(
    name='gff2neo',
    version='0.1',
    description='Parses M.tuberculosis annotation (in GFF file and online sources) '
                'and builds a Neo4j graph database storing the annotate features. '
                'It also maps these features to external services such as UniProt, CheMBL, DrugBank, etc.',
    keywords='tuberculosis, neo4j, bioservices, gff',
    packages=['gff2neo'],
    install_requires=[
        'click',
        'bioservices',
        'bcbio-gff',
        'biopython',
        'beautifulsoup4'
    ],
    entry_points={
        'console_scripts': ['gff2neo=gff2neo.cli:cli']
    },
)
