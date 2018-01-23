#!/usr/bin/env bash

GFF_FILES_DIR=data/gff_files

gff2neo examine_gff --gff_files "${GFF_FILES_DIR}"
gff2neo load_gff --gff_files "${GFF_FILES_DIR}"
gff2neo load_uniprot_data --gff_files "${GFF_FILES_DIR}"
gff2neo load_drugbank_data --gff_files "${GFF_FILES_DIR}"
gff2neo load_go_terms --gff_files "${GFF_FILES_DIR}"
gff2neo load_publications --gff_files "${GFF_FILES_DIR}"
gff2neo load_kegg_pathways --gff_files "${GFF_FILES_DIR}"
gff2neo load_reactome_pathways --gff_files "${GFF_FILES_DIR}"
