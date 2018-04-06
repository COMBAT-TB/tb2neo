from setuptools import setup

setup(
    name='gff2neo',
    version='0.0.5',
    description='Parses M.tuberculosis annotation (in GFF file and online sources) '
                'and builds a Neo4j graph database storing the annotate features. '
                'It also maps these features to external services such as UniProt, CheMBL, DrugBank, etc.',
    keywords='tuberculosis, neo4j, bioservices, gff',
    license="MIT",
    packages=['gff2neo'],
    py_modules=['gff2neo'],
    include_package_data=True,
    package_data={
        'gff2neo': ['data/drugbank/*.csv', 'data/string/*.txt', 'data/uniprot/*.csv', 'data/tbdtdb/*.txt'],
    },
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
