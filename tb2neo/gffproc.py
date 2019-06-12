"""
Interface for GFF processing.
"""
from __future__ import print_function

import os
import pprint
import sys

from BCBio import GFF
from BCBio.GFF import GFFExaminer
from tqdm import tqdm

from tb2neo.dbconn import (create_cds_nodes, create_featureloc_nodes,
                           create_gene_nodes, create_mrna_nodes,
                           create_pseudogene_nodes, create_rna_nodes,
                           create_transcript_nodes, graph, map_to_location)
from combattbmodel.core import Organism

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
MYCO_GFF = os.path.join(
    CURR_DIR, "data/myco/Mycobacterium_tuberculosis_H37Rv.gff")


def examine_gff_file(gff_file):
    """
    Examine GFF file
    :param gff_file:
    :return:
    """
    examiner = GFFExaminer()
    in_file = open(gff_file)
    pprint.pprint(examiner.available_limits(in_file))
    in_file.close()


def parse_gff(gff_file):
    """
    Parse GFF file
    :return:
    """
    # organism = create_organism_nodes(gff_file)
    organism = Organism.select(graph).first()
    # we are not interested in exons as this is a bacterial genome
    limits = [["mRNA"], ["gene", "CDS", "tRNA_gene",
                         "ncRNA_gene", "rRNA_gene", "rRNA"], ["pseudogene"]]
    for limit in limits:
        sys.stdout.write("\nLoading {}...".format(limit))
        # print("\nLoading", limit, "...")
        load_gff_data(gff_file, limit, organism)
    # print("Done.")
    sys.stdout.write("Done.")
    return gff_file


def get_locus_tags(gff_file, chunk):
    """
    Return a list of locus tags from gff_file
    :param gff_file:
    :param chunk
    :return:
    """
    sys.stdout.write("Getting locus_tags from {}...\n".format(gff_file))
    count = 0
    locus_tags = []
    for rec in GFF.parse(gff_file,
                         limit_info=dict(gff_type=['gene', 'pseudogene'])):
        for gene in rec.features:
            locus_tag = gene.qualifiers.get("gene_id", " ")[0]
            count += 1
            locus_tags.append(locus_tag)
            if count == chunk:
                yield locus_tags
                locus_tags = []
                count = 0
    yield locus_tags


def load_gff_data(gff_file, limit, organism):
    """
    Extract and load features to Neo4j
    :param organism:
    :param gff_file:
    :param limit:
    :return:
    """
    in_file = open(gff_file)
    limit_info = dict(gff_type=limit)
    for rec in GFF.parse(gff_file, limit_info=limit_info):
        for feature in tqdm(rec.features):
            rna = ["tRNA_gene", "ncRNA_gene", "rRNA_gene"]
            create_featureloc_nodes(feature)
            if feature.type == 'gene':
                create_gene_nodes(feature, organism)
                map_to_location(feature)
            elif feature.type == 'pseudogene':
                create_pseudogene_nodes(feature, organism)
                map_to_location(feature)
            elif feature.type == 'CDS':
                create_cds_nodes(feature)
                map_to_location(feature)
            elif feature.type == 'mRNA':
                create_mrna_nodes(feature)
                map_to_location(feature)
            elif feature.type in rna:
                create_rna_nodes(feature)
                map_to_location(feature)
            elif feature.type == 'transcript':
                create_transcript_nodes(feature)
                map_to_location(feature)

    in_file.close()


def map_functional_category(gff=None):
    """
    Adding functional categories to Feature Nodes
    :param gff:
    :return:
    """
    sys.stdout.write("\nAdding functional categories...")
    with open(gff) as myco_gff:
        for line in myco_gff:
            if 'Functional_Category' in str(line):
                tab_split = line.split('\t')
                # start = tab_split[3]
                end = int(tab_split[4])
                info = tab_split[8].split(";")
                functional_category = info[-1].split("=")[-1].strip()
                result = graph.run(
                    f"match(n)-[]-(l:Location) where l.fmax = {end} return n"
                ).data()
                if result:
                    for item in result:
                        node = item['n']
                        node['category'] = functional_category
                        # update RNA nodes
                        if node['biotype'] and 'rna' in str(
                                node['biotype']).lower():
                            node['description'] = info[3].replace(
                                'Product=', '')
                        node.push()

# def scrap_tbdtdb(locus_tags):
#     "http://www.bioinformatics.org/tbdtdb/tdetail.php?tid=Rv0005"
#     url = "http://www.bioinformatics.org/tbdtdb/tdetail.php?tid="
#     with open("tbdtdb_res.txt", "w") as tbdtdb:
#         for tag_list in locus_tags:
#             for lt in tqdm(tag_list):
#                 html = urlopen(url + "{}".format(lt))
#                 bs_obj = BeautifulSoup(html.read(), "html.parser")
#                 links = (bs_obj.find_all(href=re.compile('ddetail')))
#                 if len(links) > 0:
#                     for drug in links:
#                         print("\n", lt, "is targeted by", drug.text)
#                         tbdtdb.write("{}\t{}\n".format(lt, drug.text))
