from __future__ import print_function
import sys
import re
import requests
from bs4 import BeautifulSoup


def fetch_ortholog(locus_tag, ortholog_type='cdc1551'):
    ortholog_name = None
    if ortholog_type == 'cdc1551':
        tuberculist_url = \
            'http://tuberculist.epfl.ch/quicksearch.php?gene+name='
        response = requests.get(tuberculist_url + locus_tag)
        bs_obj = BeautifulSoup(response.text, 'html.parser')
        links = (bs_obj).find_all(href=re.compile('cmr.jcvi.org'))
        ortholog_name = links[0].text if links else None
    else:
        print('unknown ortholog type:', type, file=sys.stderr)
    return ortholog_name
