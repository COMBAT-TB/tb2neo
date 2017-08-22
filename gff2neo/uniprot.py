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
    :param taxonomy:
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


def write_to_csv(results):
    print("\nWriting to csv...")
    fieldnames = [
        'Entry', 'Entry_Name', 'Gene_Names_OL', 'Gene_Name', 'GO_IDs', 'InterPro', 'Interacts_With',
        'Gene_Names_Prim', 'Domain_FT', 'Protein_Names', 'GO', 'PubMed', '3D', 'Function_CC', 'Sequence',
        'Mass', 'Length', 'Protein_Families', 'GO_BP', 'GO_MF', 'GO_CC', 'Gene_SYN', 'Gene_Name_ORF', 'SeqVersion'
    ]
    print(len(fieldnames))
    with open("data/uniprot_data.csv", "w") as csv_file:
        writer = csv.DictWriter(csv_file, delimiter=',', fieldnames=fieldnames)
        writer.writeheader()
        # writer = csv.writer(csv_file)
        mapped_list = []
        for row in results:
            inner_dict = dict(zip(fieldnames, row))
            mapped_list.append(inner_dict)
        writer.writerows(mapped_list)


def query_uniprot(locus_tags, taxonomy='83332', proteome='UP000001584'):
    """
    Get data from UniProt
    :param proteome:
    :param taxonomy:
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
    # print("\nWriting to csv...")
    # with open("data/uniprot_data.csv", "w") as csv_file:
    #     writer = csv.writer(csv_file)
    for tag_list in locus_tags:
        query = '(' + '+OR+'.join(['gene:' + name for name in tag_list]) + ')'
        result = search_uniprot(query, columns, taxonomy=taxonomy,
                                proteome=proteome)
        # writer.writerows(result)
        uniprot_data.append(result)

    for data in uniprot_data:
        for entry in data:
            results.append(entry)
    end = time()
    print("\nDone fetching data from UniProt in ", end - start, "secs.")
    write_to_csv(results)
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
