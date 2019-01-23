from __future__ import print_function

import re

import requests
from bs4 import BeautifulSoup

TUBERCULIST_URL = 'http://svitsrv8.epfl.ch/tuberculist/quicksearch.php?gene+name='
STRAIN = 'CDC1551'


# Scrape CDC1551 Orthologs from Tuberculist
def fetch_ortholog(locus_tag=None, strain=STRAIN, tuberculist_url=TUBERCULIST_URL):
    """
    Orthologs
    :param locus_tag:
    :param strain:
    :param tuberculist_url:
    :return:
    """
    if locus_tag:
        try:
            response = requests.get(tuberculist_url + locus_tag)
            bs_obj = BeautifulSoup(response.text, 'html.parser')
            links = bs_obj.find_all(href=re.compile('cmr.jcvi.org'))
            ortholog_name = links[0].text if links else None
        except Exception as e:
            raise e
    else:
        raise ValueError(
            "Something is missing. locus_tag = {} and strain = {}".format(locus_tag, strain))
    return ortholog_name
