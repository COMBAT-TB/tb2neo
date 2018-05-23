#!/usr/bin/env bash

GFF_FILES_DIR=data/gff_files
OPERON_DIR=data/operons

gff2neo examine_gff "${GFF_FILES_DIR}"
# Delete existing data
gff2neo delete
gff2neo load_chromosome h37rv
gff2neo load_gff "${GFF_FILES_DIR}"
#gff2neo load_operons "${OPERON_DIR}"
gff2neo load_uniprot_data "${GFF_FILES_DIR}"
gff2neo load_drugbank_data
gff2neo load_go_terms
gff2neo load_publications
gff2neo load_kegg_pathways
gff2neo load_reactome_pathways
