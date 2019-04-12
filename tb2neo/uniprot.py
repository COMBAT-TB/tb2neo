"""
Interface to the `UniProt <http://www.uniprot.org>`_ service.
"""

from __future__ import print_function

import csv
import os
import sys
from io import StringIO

from bioservices import UniProt

uniprot_ = UniProt(verbose=False)
uniprot_._url = "https://www.uniprot.org"
CURR_DIR = os.path.dirname(os.path.abspath(__file__))

UNIPROT_DATA = os.path.join(CURR_DIR, "data/uniprot/uniprot_data.csv")


def search_uniprot(query, columns, taxonomy, proteome):
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

    result = uniprot_.search(query=query, frmt="tab",
                             columns=columns, sort=None)
    reader = csv.reader(StringIO(result), delimiter='\t')
    try:
        next(reader)
    except StopIteration:
        return []
    else:
        return list(reader)


def write_to_csv(results):
    """
    Write to CSV
    :param results:
    :return:
    """
    file_exists = os.path.isfile(UNIPROT_DATA)
    sys.stdout.write("\nWriting to csv...")
    fieldnames = [
        'Entry', 'Entry_Name', 'Gene_Names_OL', 'Gene_Name', 'GO_IDs',
        'InterPro', 'Interacts_With',
        'Gene_Names_Prim', 'Domain_FT', 'Protein_Names', 'GO', 'PubMed', '3D',
        'Function_CC', 'Sequence',
        'Mass', 'Length', 'Protein_Families', 'GO_BP', 'GO_MF', 'GO_CC',
        'Gene_SYN', 'Gene_Name_ORF',
        'SeqVersion'
    ]
    with open(UNIPROT_DATA, "a") as csv_file:
        writer = csv.DictWriter(csv_file, delimiter=',', fieldnames=fieldnames)
        if not file_exists:
            writer.writeheader()
        mapped_list = []
        for row in results:
            fn_row = dict(zip(fieldnames, row))
            mapped_list.append(fn_row)
        writer.writerows(mapped_list)
    sys.stdout.write("\nDone Writing to CSV...")


# TODO Change tax and proteome based on strain
def query_uniprot(locus_tags, taxonomy, proteome):
    """
    Get data from UniProt
    :param proteome:
    :param taxonomy:
    :param locus_tags:
    :return:
    """
    sys.stdout.write("\nQuerying UniProt...")
    uniprot_data = []
    results = []
    columns = "id, entry name, genes(OLN), genes, go-id, interpro, " \
              "interactor, genes(PREFERRED), feature(DOMAIN EXTENT), " \
              "protein names, go, citation, 3d, comment(FUNCTION), " \
              "sequence, mass, length, families, go(biological process), " \
              "go(molecular function), go(cellular component), " \
              " genes(ALTERNATIVE), genes(ORF), version(sequence)"
    # print("\nWriting to csv...")
    # with open("data/uniprot_data_0.csv", "w") as csv_file:
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
    sys.stdout.write("Done: UniProt data.")
    if len(results) > 0:
        write_to_csv(results)
    return results


def eu_mapping(from_, to):
    """
    Mapping UniProt entry to external databases
    :param to:
    :param from_:
    :return:
    """
    xref_id = None
    if from_ and to:
        _map = uniprot_.mapping(fr='ID', to=to, query=from_)
        if len(_map) != 0:
            xref_id = _map[from_]
    else:
        raise ValueError("Can't map {} to {}".format(from_, to))
    return xref_id
