"""
Interface to NCBI.
"""
import time

try:
    from urllib2 import HTTPError
except ImportError:
    from urllib.error import HTTPError

from Bio import Entrez
from Bio import Medline


# def fetch_publications(citation):
#     """
#     Fetch Publications.
#     :param citation:
#     :return:
#     """
#     print("=========================================")
#     print("Fetching publication data from PubMed.")
#     print("=========================================")
#     time.sleep(2)
#     Entrez.email = 'A.N.Other@example.com'
#     try:
#         h = Entrez.efetch(db='pubmed', id=citation, rettype='medline', retmode='text')
#     except HTTPError:
#         time.sleep(200)
#         h = Entrez.efetch(db='pubmed', id=citation, rettype='medline', retmode='text')
#     records = Medline.parse(h)
#     return records


def fetch_publication_list(citations, rettype='medline'):
    """
    Fetch Publications.
    :param rettype:
    :param citations:
    :return:
    """
    print("=====================================")
    print("Fetching {} publications from PubMed.".format(len(citations)))
    print("=====================================")
    citation_string = ','.join(citations)
    Entrez.email = 'support@sanbi.ac.za'
    retries = 5
    failed = True
    for i in range(retries):
        try:
            h = Entrez.efetch(db='pubmed', id=citation_string, rettype=rettype, retmode='text')
            failed = False
        except HTTPError:
            pass
        else:
            break
        finally:
            time.sleep(0.4)  # we are not allowed to hit NCBI more than 3 times per second
    if failed:
        print("Retrieval from PubMed failed")
        records = []
    else:
        if rettype == 'medline':
            records = Medline.parse(h)
        else:
            records = Entrez.parse(h)
    return records


def get_fasta():
    Entrez.email = 'A.N.Other@example.com'
    organism = 'Mycobacterium_tuberculosis_H37Rv'
    accession = "NC_000962.3"
    # gene_query = "{gene}[Gene] AND {organism}[Organism] AND {accession}[Accession] " \
    #              "AND RefSeq[Filter]".format(gene="", organism=organism, accession=accession)
    q = "{organism}[Organism] AND {accession}[Accession] AND RefSeq[Filter]".format(organism=organism,
                                                                                    accession=accession)
    handle = Entrez.esearch(db="nucleotide", term=q)
    record = Entrez.read(handle)
    ids = record[u'IdList']
    seq_id = ids[0]  # you must implement an if to deal with <0 or >1 cases
    handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
    record = handle.read()
    # remove >NC_000962.3 Mycobacterium tuberculosis H37Rv, complete genome
    residues = ''.join(record.strip('\n').split('\n')[1:])
    return residues
