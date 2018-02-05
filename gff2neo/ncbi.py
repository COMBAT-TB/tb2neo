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
