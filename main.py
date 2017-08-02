#!/usr/bin/env python
from gff2neo.gffproc import *
from gff2neo.uniprot import *

gff_file = "data/MTB_H37rv.gff3"

if __name__ == '__main__':
    import time

    time.sleep(10)
    delete_data()
    examine(gff_file)
    parse_gff(gff_file)
    # build_relationships()
    # query_uniprot(get_locus_tags(gff_file, 400))
