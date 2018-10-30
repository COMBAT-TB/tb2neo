"""
Interface to NCBI.
"""
import sys
import time

try:
    from urllib2 import HTTPError
except ImportError:
    from urllib.error import HTTPError

from Bio import Entrez
from Bio import Medline


def search_pubmed(genename):
    """
    Search PubMed
    :param genename:
    :return:
    """
    term = "{genename} AND tuberculosis".format(genename=genename)
    Entrez.email = 'support@sanbi.ac.za'
    try:
        h = Entrez.esearch(db='pubmed', term=term)
        result = Entrez.read(h)
        ids = result['IdList']
    except HTTPError as e:
        sys.stderr.write("{e}".format(e=e))
    return ids


def fetch_publication_list(citations, rettype='medline'):
    """
    Fetch Publications.
    :param rettype:
    :param citations:
    :return:
    """
    sys.stdout.write("=====================================")
    print("Fetching {} publications. rettype: {}."
          .format(len(citations), rettype))
    citation_string = ','.join(citations)
    Entrez.email = 'support@sanbi.ac.za'
    retries = 5
    failed = True
    for i in range(retries):
        try:
            h = Entrez.efetch(db='pubmed', id=citation_string,
                              rettype=rettype, retmode='text')
            failed = False
        except HTTPError:
            pass
        else:
            break
        finally:
            # we are not allowed to hit NCBI more than 3 times per second
            time.sleep(0.4)
    if failed:
        print("Retrieval from PubMed failed")
        records = []
    else:
        if rettype == 'medline':
            records = Medline.parse(h)
        else:
            records = Entrez.parse(h)
    return records


def get_fasta(strain):
    Entrez.email = 'support@sanbi.ac.za'
    if strain and strain == "h37rv":
        organism = 'Mycobacterium_tuberculosis_H37Rv'
        accession = "NC_000962.3"
    elif strain and strain == "cdc1551":
        organism = 'Mycobacterium_tuberculosis_CDC1551'
        accession = "NC_002755.2"
    else:
        organism = None
        accession = None
    # gene_query = "{gene}[Gene] AND {organism}[Organism] " \
    #              "AND {accession}[Accession] AND RefSeq[Filter]" \
    #              "".format(gene="", organism=organism, accession=accession)
    q = "{organism}[Organism] AND {accession}[Accession] " \
        "AND RefSeq[Filter]".format(organism=organism,
                                    accession=accession)
    handle = Entrez.esearch(db="nucleotide", term=q)
    record = Entrez.read(handle)
    ids = record[u'IdList']
    seq_id = ids[0]  # you must implement an if to deal with <0 or >1 cases
    handle = Entrez.efetch(db="nucleotide", id=seq_id,
                           rettype="fasta", retmode="text")
    record = handle.read()
    # remove >NC_000962.3 Mycobacterium tuberculosis H37Rv, complete genome
    residues = ''.join(record.strip('\n').split('\n')[1:])
    if isinstance(residues, str):
        return residues
    else:
        return ''
