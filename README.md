# tb2neo

[![Build Status](https://travis-ci.org/COMBAT-TB/tb2neo.svg?branch=master)](https://travis-ci.org/COMBAT-TB/tb2neo) [![Coverage Status](https://coveralls.io/repos/github/COMBAT-TB/tb2neo/badge.svg?branch=master)](https://coveralls.io/github/COMBAT-TB/tb2neo?branch=master)

Parses _M.tuberculosis_ annotation (GFF file) and builds a Neo4j graph
database storing the annotated features. _tb2neo_ also aggregates and maps
these annotated features to external services such as UniProt, CheMBL,
DrugBank, KEGG, Reactome, QuickGo, STRING-DB etc.

## Usage

### Prerequisite

- [Docker](https://docs.docker.com/) or [Neo4j installation](https://neo4j.com/docs/getting-started/current/get-started-with-neo4j/)
- `python-pip`

**Pull and run the [neo4j docker image](https://hub.docker.com/_/neo4j/):**

```sh
$ docker run -d -p 7474:7474 -p 7687:7687 --name neo -e NEO4J_AUTH=none -v=$HOME/neo4j/data:/data neo4j:3.5
```

**Clone repository and create a virtual environment:**

```sh
$ git clone https://github.com/COMBAT-TB/tb2neo.git
$ cd tb2neo
$ virtualenv envname
$ source envname/bin/activate
```

**Install and run `tb2neo`:**

```sh
$ pip install -r requirements.txt
$ python setup.py install
$ tb2neo --help
$ tb2neo load_gff --gff_files PATH/TO/TB_GFF3_FILES
```

Point your browser to [localhost:7474](http://localhost:7474]).

### `db.schema()`

![DB_MODEL](https://raw.githubusercontent.com/COMBAT-TB/tb2neo/master/images/dbschema_.png)
