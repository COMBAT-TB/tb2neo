from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='tb2neo',
    version='0.0.6',
    url='https://github.com/COMBAT-TB/tb2neo',
    description='Builds a M.tb annotation graph database from GFF files',
    long_description=long_description,
    long_description_content_type="text/markdown",
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
