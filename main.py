#!/usr/bin/env python
from gff2neo.gffproc import *

gff_file = "data/MTB_H37rv.gff3"


def main():
    if os.path.isdir(os.getcwd() + "/data") and os.path.exists(gff_file):
        import time
        time.sleep(10)
        delete_data()
        examine(gff_file)
        parse_gff(gff_file)
        build_relationships()
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
