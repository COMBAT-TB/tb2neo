#!/usr/bin/env python
from gff2neo.gffproc import *

gff_file = "data/MTB_H37rv.gff3"


def main():
    if os.path.isdir(os.getcwd() + "/data") and os.path.exists(gff_file):
        import time
        time.sleep(10)
        delete_data()
        examine(gff_file)
        gff_p_st = time.time()
        parse_gff(gff_file)
        gff_p_end = time.time()
        sys.stdout.write("Done parsing GFF in {} seconds.".format(gff_p_end - gff_p_st))
        rel_st = time.time()
        build_relationships()
        rel_end = time.time()
        sys.stdout.write("Built GFF relationships in {} seconds.".format(rel_end - rel_st))
        if os.path.exists(uniprot_data_csv) and os.stat(uniprot_data_csv).st_size > 0:
            create_uniprot_nodes()
            create_go_term_nodes()
        else:
            query_uniprot(get_locus_tags(gff_file, 400))
            create_uniprot_nodes()
            create_go_term_nodes()
    else:
        raise Exception("Couldn't find H37Rv GFF file: {}!".format(gff_file))


if __name__ == '__main__':
    main()
