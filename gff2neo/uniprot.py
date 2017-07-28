"""
Interface to the `UniProt <http://www.uniprot.org>`_ service.
"""

from __future__ import print_function
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import csv
from time import time
from bioservices import UniProt

u = UniProt(verbose=False)


def search_uniprot(query, columns, taxonomy='83332', proteome='UP000001584'):
    """
    Search UniProt and return results as list
    :param query:
    :param columns:
    :param proteome:
    :return:
    """
    query = "taxonomy:{}+AND+proteome:{}+AND+{}".format(taxonomy,
                                                        proteome, query)

    result = u.search(query=query, frmt="tab", columns=columns, sort=None)
    reader = csv.reader(StringIO(result), delimiter='\t')
    try:
        next(reader)
    except StopIteration:
        return []
    else:
        return list(reader)


def query_uniprot(locus_tags, taxonomy='83332', proteome='UP000001584'):
    """
    Get data from UniProt
    :param locus_tags:
    :return:
    """
    print("Querying UniProt...")
    start = time()
    columns = "id, entry name, genes(OLN), genes, go-id, interpro, " \
              "interactor, genes(PREFERRED), feature(DOMAIN EXTENT), " \
              "protein names, go, citation, 3d, comment(FUNCTION), " \
              "sequence, mass, length, families, go(biological process), " \
              "go(molecular function), go(cellular component), " \
              " genes(ALTERNATIVE), genes(ORF), version(sequence)"
    uniprot_data = []
    results = []
    for tag_list in locus_tags:
        query = '(' + '+OR+'.join(['gene:' + name for name in tag_list]) + ')'
        result = search_uniprot(query, columns, taxonomy=taxonomy,
                                proteome=proteome)
        for row in result:
            print(row[9], row[0], row[2], row[3])
        uniprot_data.append(result)

    for data in uniprot_data:
        for entry in data:
            results.append(entry)
    end = time()
    print("Done fetching data from UniProt in ", end - start, "secs.")
    return results


def map_ue_to_pdb(ue):
    """
    Mapping UniProt entry to PDB
    :param ue:
    :return:
    """
    pdb_id = None
    _pdb = u.mapping(fr='ID', to='PDB_ID', query=ue)
    if len(_pdb) != 0:
        pdb_id = _pdb[ue]
    return pdb_id
