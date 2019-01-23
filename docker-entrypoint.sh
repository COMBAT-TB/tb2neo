#!/usr/bin/env bash

GFF_FILES_DIR=data/gff_files
OPERON_DIR=data/operons
MUTATIONS_DIR=tb2neo/data/mutations

#tb2neo examine_gff "${GFF_FILES_DIR}"
# Delete existing data
tb2neo delete
tb2neo load_organism "${GFF_FILES_DIR}"
tb2neo load_chromosome h37rv
tb2neo load_gff "${GFF_FILES_DIR}"
# tb2neo load_operons "${OPERON_DIR}"
tb2neo load_uniprot_data "${GFF_FILES_DIR}"
tb2neo load_drugbank_data
tb2neo load_known_mutations "${MUTATIONS_DIR}"
tb2neo load_go_terms
tb2neo load_kegg_pathways
tb2neo load_reactome_pathways
tb2neo load_publications

