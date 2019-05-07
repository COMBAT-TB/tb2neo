from setuptools import find_packages, setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='tb2neo',
    version='0.0.8',
    url='https://github.com/COMBAT-TB/tb2neo',
    description='Builds a M.tb annotation graph database from GFF files',
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords='tuberculosis, neo4j, bioservices, gff',
    license="GPLv3",
    py_modules=['tb2neo'],
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'tb2neo': ['data/drugbank/*.csv', 'data/string/*.txt',
                   'data/uniprot/*.csv', 'data/tbdtdb/*.txt'],
    },
    python_requires='~=3.6',
    install_requires=[
        'click',
        'py2neo==3.1.2',
        'bioservices',
        'bcbio-gff',
        'biopython',
        'beautifulsoup4'
    ],
    entry_points={
        'console_scripts': ['tb2neo=tb2neo.cli:cli']
    },
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Lavnguage :: Python',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
)
