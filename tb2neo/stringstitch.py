"""Fetching data from STRING and STITCH
"""
import requests
from pandas import read_csv
from requests.exceptions import ConnectionError, HTTPError

from tb2neo.dbconn import graph
from tb2neo.uniprot import UNIPROT_DATA, eu_mapping

STRING_API_URL = "https://string-db.org/api"
STITCH_API_URL = "http://stitch.embl.de/api"
TAXON_ID = "83332"


def create_ppi(rv_a, rv_b, score):
    cypher_q = f"MATCH (ga)--(ap:Protein),(gb)--(bp:Protein) "
    cypher_q += f"WHERE ga.uniquename = '{rv_a}' AND gb.uniquename = '{rv_b}' "
    cypher_q += f"CREATE (ap)-[r:INTERACTS_WITH {{ score: {score} }}]->(bp)"
    print(cypher_q)
    data = graph.run(cypher_q).data()
    return data


def fetch_string_data(gene, output_format='json', method='network',
                      taxon=TAXON_ID):
    request_url = STRING_API_URL + "/" + output_format + "/" + method + "?"
    request_url += f"identifiers={gene}"
    request_url += "&" + f"species={TAXON_ID}"
    result = None
    try:
        response = requests.get(request_url)
    except (HTTPError, ConnectionError) as error:
        print(f'An error occured: {error}')
    except Exception as e:
        print(f'{e}')
    else:
        if response.status_code == 200:
            print(request_url)
            result = response.json()
        return result


def load_string_data():
    df = read_csv(UNIPROT_DATA).fillna("")
    for entry in df.values:
        gene = eu_mapping(from_=entry[0], to='TUBERCULIST_ID')[0]
        # This is interesting
        # https://string-db.org/api/json/network?identifiers=None&species=83332
        data = fetch_string_data(gene=gene)
        for ppi in data:
            create_ppi(ppi["stringId_A"], ppi["stringId_B"], ppi["score"])


def fetch_stitch_data():
    pass
