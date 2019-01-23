from setuptools import setup

setup(
    name='tb2neo',
    version='0.0.5',
    description='Parses M.tuberculosis annotation (in GFF file and online sources) '
                'and builds a Neo4j graph database storing the annotate features. '
                'It also maps these features to external services such as UniProt, CheMBL, DrugBank, etc.',
    keywords='tuberculosis, neo4j, bioservices, gff',
    license="MIT",
    packages=['tb2neo'],
    py_modules=['tb2neo'],
    include_package_data=True,
    package_data={
        'tb2neo': ['data/drugbank/*.csv', 'data/string/*.txt',
                   'data/uniprot/*.csv', 'data/tbdtdb/*.txt'],
    },
    install_requires=[
        'click',
        'bioservices',
        'bcbio-gff',
        'biopython',
        'beautifulsoup4'
    ],
    entry_points={
        'console_scripts': ['tb2neo=tb2neo.cli:cli']
    },
)
