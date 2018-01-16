#!/usr/bin/env python
from gff2neo.gffproc import *

gff_file = "data/MTB_H37rv.gff3"


def check_csv(csvfile):
    """
    Check if csv file exists and is not empty
    :param csvfile:
    :return:
    """
    _file = False
    if os.path.exists(csvfile) and os.stat(csvfile).st_size > 0:
        _file = True
    return _file


def main():
    if os.path.isdir(os.getcwd() + "/data") and os.path.exists(gff_file):
        import time
        time.sleep(20)
        delete_db_data()
        examine_gff_file(gff_file)
        start_time = time.time()
        parse_gff(gff_file)
        end_time = time.time()
        sys.stdout.write("\nDone parsing GFF in {} seconds.\n".format(end_time - start_time))
        start_time = time.time()
        build_gff_relationships()
        end_time = time.time()
        sys.stdout.write("\nBuilt GFF relationships in {} seconds.\n".format(end_time - start_time))
        if check_csv(uniprot_data_csv):
            create_protein_nodes()
            create_kegg_pathways_nodes()
            create_reactome_pathway_nodes()
            create_go_term_nodes()
            if check_csv(target_protein_ids_csv) and check_csv(drug_vocab_csv):
                create_drugbank_nodes()
            create_publication_nodes()
            map_gene_to_protein(get_locus_tags(gff_file, 400))
            # build_string_ppis()
        else:
            query_uniprot(get_locus_tags(gff_file, 400))
            create_protein_nodes()
            create_kegg_pathways_nodes()
            create_reactome_pathway_nodes()
            create_go_term_nodes()
            if check_csv(target_protein_ids_csv) and check_csv(drug_vocab_csv):
                create_drugbank_nodes()
            create_publication_nodes()
            map_gene_to_protein(get_locus_tags(gff_file, 400))
            # build_string_ppis()
    else:
        raise Exception("Couldn't find H37Rv GFF file: {}!".format(gff_file))


if __name__ == '__main__':
    main()
