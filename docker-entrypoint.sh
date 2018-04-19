#!/usr/bin/env bash

GFF_FILES_DIR=data/gff_files

#gff2neo examine_gff "${GFF_FILES_DIR}"
#gff2neo load_gff "${GFF_FILES_DIR}"
#gff2neo load_uniprot_data "${GFF_FILES_DIR}"
#gff2neo load_drugbank_data
gff2neo load_go_terms
#gff2neo load_publications
#gff2neo load_kegg_pathways
#gff2neo load_reactome_pathways
