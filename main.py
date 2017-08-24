#!/usr/bin/env python
import os

from gff2neo.gffproc import *

gff_file = "data/MTB_H37rv.gff3"

if __name__ == '__main__':
    import time

    time.sleep(10)
    delete_data()
    examine(gff_file)
    parse_gff(gff_file)
    build_relationships()
    if os.path.exists(uniprot_data_csv) and os.stat(uniprot_data_csv).st_size > 0:
        create_uniprot_nodes()
    else:
        query_uniprot(get_locus_tags(gff_file, 400))
        create_uniprot_nodes()
