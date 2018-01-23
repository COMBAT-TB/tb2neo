from __future__ import print_function

import re

import requests
from bs4 import BeautifulSoup


def fetch_ortholog(locus_tag=None, strain='cdc1551'):
    # ortholog_name = None
    if locus_tag and strain == 'cdc1551':
        tuberculist_url = \
            'http://tuberculist.epfl.ch/quicksearch.php?gene+name='
        try:
            response = requests.get(tuberculist_url + locus_tag)
            bs_obj = BeautifulSoup(response.text, 'html.parser')
            links = bs_obj.find_all(href=re.compile('cmr.jcvi.org'))
            ortholog_name = links[0].text if links else None
        except Exception as e:
            raise e
    else:
        raise ValueError(
            "Something is missing. locus_tag = {} and strand = {}".format(locus_tag, strain))
    return ortholog_name
