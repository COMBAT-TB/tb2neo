#!/usr/bin/env bash

GFF_FILES_DIR=data/gff_files
OPERON_DIR=data/operons
MUTATIONS_DIR=tb2neo/data/mutations

#tb2neo examine_gff "${GFF_FILES_DIR}"
# Delete existing data
tb2neo delete
tb2neo load-organism "${GFF_FILES_DIR}"
tb2neo load-chromosome h37rv
tb2neo load-gff "${GFF_FILES_DIR}"
# tb2neo load-operons "${OPERON_DIR}"
tb2neo load-uniprot-data "${GFF_FILES_DIR}"
tb2neo load-drugbank-data
tb2neo load-known-mutations "${MUTATIONS_DIR}"
tb2neo load-go-terms
tb2neo load-kegg-pathways
tb2neo load-reactome-pathways
# tb2neo load-publications

