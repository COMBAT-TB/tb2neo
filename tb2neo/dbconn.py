"""
Interface to the Neo4j Database
"""
import csv
import os
import sys
import zipfile

from bioservices import KEGG, ChEMBL, QuickGO, ReactomeOld
from pandas import read_csv
from py2neo import Graph
from tqdm import tqdm

from tb2neo.ftpconn import get_nucleotides
from tb2neo.model.core import (CDS, Author, Chromosome, DbXref, Drug, Gene,
                               GOTerm, InterProTerm, Location, MRna, NCRna,
                               Operon, Organism, Pathway, Protein, PseudoGene,
                               Publication, RRna, Transcript, TRna)
from tb2neo.model.vcfmodel import CallSet, Variant, VariantSet
from tb2neo.ncbi import fetch_publication_list
from tb2neo.orthologs import fetch_ortholog
from tb2neo.quickgo import query_quickgo
from tb2neo.uniprot import UNIPROT_DATA, eu_mapping

graph = Graph(host=os.environ.get("DATABASE_URL", "localhost"), bolt=True,
              password=os.environ.get("NEO4J_PASSWORD", ""))

chembl = ChEMBL(verbose=False)
quick_go = QuickGO(verbose=False)
quick_go.url = 'http://www.ebi.ac.uk/QuickGO-Old'
reactome_old = ReactomeOld(verbose=False)
kegg = KEGG(verbose=False)

# watch("neo4j.bolt")

gene_dict = dict()
mrna_dict = dict()
transcript_dict = dict()
pseudogene_dict = dict()
cds_dict = dict()
exon_dict = dict()
rrna_dict = dict()
trna_dict = dict()
ncrna_dict = dict()
location_dict = dict()
go_term_set = set()

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
TARGET_PROTEIN_IDS = os.path.join(
    CURR_DIR, "data/drugbank/drugbank_approved_target_polypeptide_ids.csv.zip")
DRUG_VOCAB = os.path.join(
    CURR_DIR, "data/drugbank/drugbank_all_drugbank_vocabulary.csv.zip")
STRING_DATA = os.path.join(
    CURR_DIR, "data/string/83332.protein.links.detailed.v10.5.txt")


def delete_db_data():
    """
    Delete existing data.
    :return:
    """
    sys.stdout.write(
        "Deleting all nodes and relationships in {}\n".format(graph))

    graph.delete_all()


def create_organism_nodes(gff_file=None):
    """
    Create Organism Nodes
    :return:
    """
    # TODO: Change the strain to accommodate MTB strains
    # get strain name from gff_file name
    gff_strain = None
    if "/" in str(gff_file):
        gff_strain = str(gff_file).split("/")[-1]
    if "." in gff_strain:
        gff_strain = gff_strain.split(".")[0]
    else:
        gff_strain = gff_strain

    strain = gff_strain
    genus = "Mycobacterium"
    species = "M. tuberculosis"
    common_name = "TB"

    organism = Organism(strain=strain, genus=genus,
                        species=species, common_name=common_name)
    graph.create(organism)
    return organism


def create_chromosome_nodes(strain):
    """
    Create Chromosome Nodes
    :return:
    """
    name = "Chr1"
    uniquename = "Chr1"
    chromosome = Chromosome()
    chromosome.name = name
    chromosome.uniquename = uniquename
    chromosome.residues = get_nucleotides(strain=strain)
    graph.create(chromosome)
    organism = Organism.select(graph).first()
    chromosome.belongs_to.add(organism)
    graph.push(chromosome)


def create_gene_nodes(feature, organism):
    """
    Create Gene Nodes
    :param organism:
    :param feature:
    :return:
    """
    names = get_feature_name(feature)
    name = names.get("Name", names.get("UniqueName"))
    unique_name = names.get("UniqueName", name)
    description = feature.qualifiers.get("description", "")
    biotype = feature.qualifiers['biotype'][0]
    parent = feature.qualifiers.get("Parent", " ")[0]

    gene = Gene()
    gene.name = name
    gene.uniquename = unique_name
    gene.parent = parent[parent.find(':') + 1:]
    gene.biotype = biotype
    gene.description = description
    graph.create(gene)
    gene.belongs_to.add(organism)
    graph.push(gene)
    gene_dict[unique_name] = gene


def create_mrna_nodes(feature):
    """
    Create mRNA Nodes
    :param feature:
    :return:
    """
    names = get_feature_name(feature)
    name = names.get("Name", names.get("UniqueName"))
    unique_name = names.get("UniqueName", name)
    biotype = feature.qualifiers['biotype'][0]
    parent = feature.qualifiers.get("Parent", " ")[0]

    mrna = MRna()
    mrna.name = name
    mrna.parent = parent[parent.find(':') + 1:]
    mrna.uniquename = unique_name
    mrna.biotype = biotype
    graph.create(mrna)
    mrna_dict[unique_name] = mrna


def create_transcript_nodes(feature):
    """
    Create Transcript Nodes
    :param feature:
    :return:
    """
    names = get_feature_name(feature)
    name = names.get("Name", names.get("UniqueName"))
    unique_name = names.get("UniqueName", name)
    biotype = feature.qualifiers['biotype'][0]
    parent = feature.qualifiers.get("Parent", " ")[0]

    transcript = Transcript()
    transcript.name = name
    transcript.parent = parent[parent.find(':') + 1:]
    transcript.uniquename = unique_name
    transcript.biotype = biotype
    graph.create(transcript)
    # transcript_dict[unique_name] = transcript


def create_pseudogene_nodes(feature, organism):
    """
    Create Pseudogene Nodes
    :param organism:
    :param feature:
    :return:
    """
    names = get_feature_name(feature)
    name = names.get("Name", names.get("UniqueName"))
    unique_name = names.get("UniqueName", name)
    description = feature.qualifiers.get("description", " ")[0]
    biotype = feature.qualifiers['biotype'][0]
    parent = feature.qualifiers.get("Parent", " ")[0]

    pseudogene = PseudoGene()
    pseudogene.name = name
    pseudogene.uniquename = unique_name
    pseudogene.parent = parent[parent.find(':') + 1:]
    pseudogene.description = description
    pseudogene.biotype = biotype
    graph.create(pseudogene)
    pseudogene.belongs_to.add(organism)
    graph.push(pseudogene)
    pseudogene_dict[unique_name] = pseudogene


def create_rna_nodes(feature):
    """
    Create RNA Nodes
    :param feature:
    :return:
    """
    names = get_feature_name(feature)
    name = names.get("Name", names.get("UniqueName"))
    unique_name = names.get("UniqueName", name)
    biotype = feature.qualifiers['biotype'][0]
    parent = feature.qualifiers.get("Parent", " ")[0]

    if feature.type == 'tRNA_gene' and biotype == "tRNA":
        trna = TRna()
        trna.name = name
        trna.parent = parent[parent.find(':') + 1:]
        trna.uniquename = unique_name
        trna.biotype = biotype
        graph.create(trna)
        trna_dict[unique_name] = trna
    if feature.type == 'ncRNA_gene' and biotype == "ncRNA":
        ncrna = NCRna()
        ncrna.name = name
        ncrna.parent = parent[parent.find(':') + 1:]
        ncrna.uniquename = unique_name
        ncrna.biotype = biotype
        graph.create(ncrna)
        ncrna_dict[unique_name] = ncrna
    if feature.type == 'rRNA_gene' and biotype == "rRNA":
        rrna = RRna()
        rrna.name = name
        rrna.parent = parent[parent.find(':') + 1:]
        rrna.uniquename = unique_name
        rrna.biotype = biotype
        graph.create(rrna)
        rrna_dict[unique_name] = rrna


def create_cds_nodes(feature):
    """
    Create CDS Nodes
    :param feature:
    :return:
    """
    names = get_feature_name(feature)
    name = names.get("Name", names.get("UniqueName"))
    unique_name = names.get("UniqueName", name)
    parent = feature.qualifiers.get("Parent", " ")[0]

    cds = CDS()
    cds.name = name
    cds.parent = parent[parent.find(':') + 1:]
    cds.uniquename = unique_name
    graph.create(cds)
    cds_dict[unique_name] = cds


def create_featureloc_nodes(feature):
    """
    Create FeatureLoc Nodes
    :param feature:
    :return:
    """
    srcfeature_id = get_feature_name(feature).get("UniqueName")
    # Add 1 to start. Ensembl GFF in one based.
    primary_key = (feature.location.start + 1) + feature.location.end
    feature_loc = Location(pk=primary_key, fmin=(feature.location.start + 1),
                           fmax=feature.location.end,
                           strand=feature.location.strand)
    graph.create(feature_loc)
    location_dict[srcfeature_id] = feature_loc


def get_feature_name(feature):
    """
    Get Feature Name and UniqueName
    :param feature:
    :return:
    """
    names = dict()
    u_name = feature.id[feature.id.find(":") + 1:]
    if feature.qualifiers.get("Name"):
        name = feature.qualifiers["Name"][0]
        names["Name"] = name[name.find(':') + 1:]
        names["UniqueName"] = u_name
    else:
        names["Name"] = names["UniqueName"] = u_name
    return names


def build_gff_relationships():
    """
    Build GFF Feature relationships
    :return:
    """
    sys.stdout.write("\nBuilding GFF Relationships...\n")
    for m, mrna in mrna_dict.items():
        if mrna.parent in gene_dict.keys():
            gene = gene_dict.get(mrna.parent)
            mrna.part_of_g.add(gene)
            graph.push(mrna)
            gene.part_of.add(mrna)
            graph.push(gene)
        for p, pseudogene in pseudogene_dict.items():
            if mrna.parent == pseudogene.uniquename:
                mrna.part_of_pg.add(pseudogene)
                graph.push(mrna)
                pseudogene.part_of.add(mrna)
                graph.push(pseudogene)
        for c, cds in cds_dict.items():
            if mrna.uniquename == cds.parent:
                cds.part_of.add(mrna)
                graph.push(cds)
                mrna.part_of_cds.add(cds)
                graph.push(mrna)
    # for t, transcript in transcript_dict.items():
    #     if transcript.parent in ncrna_dict.keys():
    #         ncrna = ncrna_dict.get(transcript.parent)
    #         ncrna.part_of.add(transcript)
    #         graph.push(ncrna)
    #     if transcript.parent in trna_dict.keys():
    #         trna = trna_dict.get(transcript.parent)
    #         trna.part_of.add(transcript)
    #         graph.push(trna)
    sys.stdout.write("\nDone Building GFF Relationships...\n")


def map_to_location(feature):
    """
    Map GFF Features to Location and Organism
    :param feature:
    :return:
    """
    # Find feature location with a srcfeature_id attr.
    # matching this features uniquename and link them via
    # LOCATED_AT
    srcfeature_id = get_feature_name(feature).get("UniqueName")
    location = location_dict.get(srcfeature_id)
    chromosome = Chromosome.select(graph).first()
    rna = ["tRNA_gene", "ncRNA_gene", "rRNA_gene"]
    if location:
        if feature.type == 'gene':
            _feature = gene_dict.get(srcfeature_id)
            _feature.location.add(location)
            _feature.located_on.add(chromosome)
            _feature.residues = _feature.get_residues()
            graph.push(_feature)
        elif feature.type == 'pseudogene':
            _feature = pseudogene_dict.get(srcfeature_id)
            _feature.location.add(location)
            _feature.located_on.add(chromosome)
            _feature.residues = _feature.get_residues()
            graph.push(_feature)
        elif feature.type in rna:
            if feature.type == 'tRNA_gene':
                _feature = trna_dict.get(srcfeature_id)
                _feature.location.add(location)
                _feature.located_on.add(chromosome)
                _feature.residues = _feature.get_residues()
                graph.push(_feature)
            if feature.type == 'ncRNA_gene':
                _feature = ncrna_dict.get(srcfeature_id)
                _feature.location.add(location)
                _feature.located_on.add(chromosome)
                _feature.residues = _feature.get_residues()
                graph.push(_feature)
            if feature.type == 'rRNA_gene':
                _feature = rrna_dict.get(srcfeature_id)
                _feature.location.add(location)
                _feature.located_on.add(chromosome)
                _feature.residues = _feature.get_residues()
                graph.push(_feature)
        elif feature.type == 'CDS':
            _feature = cds_dict.get(srcfeature_id)
            _feature.location.add(location)
            _feature.located_on.add(chromosome)
            _feature.residues = _feature.get_residues()
            graph.push(_feature)
        elif feature.type == 'mRNA':
            _feature = mrna_dict.get(srcfeature_id)
            _feature.location.add(location)
            _feature.located_on.add(chromosome)
            _feature.residues = _feature.get_residues()
            graph.push(_feature)


def create_go_term_nodes():
    """
    Create GOTerm Nodes and build Protein relationships.
    :return:
    """
    sys.stdout.write("\nCreating GoTerm Nodes...\n")
    quick_go_data_dict = dict()

    with open(UNIPROT_DATA, 'r') as csv_file:
        reader = csv.DictReader(csv_file, delimiter=',')
        for entry in reader:
            protein_entry = entry['Entry']
            protein = None
            if protein_entry:
                protein = Protein.select(graph, protein_entry).first()
            go_ids = [g for g in entry['GO_IDs'].split("; ") if
                      g and g.startswith("GO:") and g != 'GO_IDs']
            go_term_set.add(go_id for go_id in go_ids)
            terms = ','.join(go_ids)
            response = query_quickgo(terms)
            if response.status_code == 200:
                data = response.json()
                for result in data['results']:
                    accession = result['id']
                    # quick_go_data[accession] = result
                    name = result['name']
                    definition = result['definition']['text']
                    ontology = result['aspect']
                    go_term = GOTerm(accession=accession, name=name.strip(),
                                     definition=definition.strip(),
                                     ontology=ontology)
                    graph.create(go_term)
                    quick_go_data_dict[accession] = result
                    if protein is not None:
                        protein.assoc_goterm.add(go_term)
                        graph.push(protein)
                        go_term.protein.add(protein)
                        graph.push(go_term)
            else:
                sys.stderr.write(
                    '\nA status of {code} occurred for {go}\n'.format(
                        code=response.status_code, go=terms))

    sys.stdout.write("\nMapping GoTerm Relations...\n")
    for term_accession, v in quick_go_data_dict.items():
        term = GOTerm.select(graph, term_accession).first()
        if term:
            if v.get('children'):
                for child in v['children']:
                    child_id = child['id']
                    child_relation = child['relation']
                    child_go = GOTerm.select(graph, child_id).first()
                    if child_go:
                        if 'is_a' in child_relation:
                            term.is_a.add(child_go)
                        elif 'part_of' in child_relation:
                            term.part_of_go.add(child_go)
                        elif 'regulates' in child_relation:
                            term.regulates.add(child_go)
                        elif 'capable_of' in child_relation:
                            term.capable_of.add(child_go)
                    graph.push(term)
    sys.stdout.write("Created {} GoTerms.".format(len(go_term_set)))


def create_interpro_term_nodes(protein, interpro_ids):
    """
    Create InterPro Term Nodes.
    :param protein:
    :param interpro_ids:
    :return:
    """
    import time
    # http://generic-model-organism-system-database.450254.n5.nabble.com/Re-GMOD-devel-Storing-Interpro-domains-in-Chado-td459778.html
    terms = [interpro_id for interpro_id in interpro_ids.split(
        "; ") if interpro_id]
    for interpro in terms:
        dbxref = DbXref(db="InterPro", accession=interpro, version=time.time())
        graph.create(dbxref)
        protein.dbxref.add(dbxref)
        interproterm = InterProTerm(accession=interpro)
        graph.create(interproterm)
        protein.assoc_intterm.add(interproterm)
        graph.push(protein)
        interproterm.assoc_protein.add(protein)
        graph.push(interproterm)


def create_author_nodes(publication, full_author):
    """
    Create Author Nodes.
    :param publication:
    :param full_author:
    :return:
    """
    # TODO: Get more info about Authors
    if full_author:
        for au in full_author:
            _author = Author()
            _author.givennames = au
            graph.create(_author)
            _author.wrote.add(publication)
            publication.author.add(_author)
            graph.push(_author)
            graph.push(publication)


def create_publication_nodes(uniprot_data):
    """
    Create Publication Nodes
    :return:
    """
    pmid_set = set()
    sys.stdout.write("\nCreating Publication Nodes...\n")
    df = read_csv(uniprot_data).fillna("")

    for entry in tqdm(df.values):
        # TODO: Search for more publications
        # locus_tag = entry[2].strip()
        # gene_name_prim = entry[7].strip()
        # genename = gene_name_prim if gene_name_prim is not '' else locus_tag
        # pmids = search_pubmed(genename)  # list
        # if pmids:
        #     pmid_set.update(pmids)
        protein_entry = str(entry[0]).strip()
        protein = None
        if protein_entry:
            protein = Protein.select(graph, protein_entry).first()
            uniprot_pmids = [p.strip() for p in entry[11].split("; ") if p]
            pmid_set.update(uniprot_pmids)
            for p_id in uniprot_pmids:
                pub = Publication()
                pub.pmid = p_id
                graph.create(pub)
                if protein:
                    protein.published_in.add(pub)
                    graph.push(protein)
    # Let's build a set before calling API and DB
    num_ids = len(pmid_set)
    sys.stdout.write(f"\nFetching data for {num_ids} publications\n")
    chunksize = 500
    records = []
    for start in range(0, num_ids, chunksize):
        subset = list(pmid_set)[start: start + chunksize]
        # Remove pmid from subset if it exists in records.
        [subset.pop(k) for k, v in enumerate(subset)
            if [rec['PMID'] for rec in records if v == rec["PMID"]]]
        records.extend(fetch_publication_list(subset))
    for record in tqdm(records):
        if isinstance(record, dict):
            # https://www.nlm.nih.gov/bsd/mms/medlineelements.html
            pm_id = record['PMID']
            # there is internal caching so using a dictionary here doesn't
            # actually seem to save any time - pvh
            title = record.get('TI')
            volume = record.get('VI')
            issue = record.get('IP')
            pages = record.get('PG')
            date_of_pub = record.get('DP')
            pub_place = record.get('PL')
            publisher = record.get('SO')
            # author = record.get('AU')
            full_author = record.get('FAU')
        else:
            pm_id = record['id:'][0][record['id:'][0].find('able: ') + 6:]
            record = fetch_publication_list(pm_id, rettype='xml')
            rec = next(record)
            article = rec['MedlineCitation']['Article']
            title = article['ArticleTitle']
            pages = article['Pagination']['MedlinePgn']
            volume = article['Journal']['JournalIssue']['Volume']
            issue = article['Journal']['JournalIssue']['Issue']
            date_of_pub = article['Journal']['JournalIssue']['PubDate']['Month'] + \
                " " + article['Journal']['JournalIssue']['PubDate']['Year']
            pub_place = rec['MedlineCitation']['MedlineJournalInfo']['Country']
            publisher = None
            # author = None
            # full_author = article['AuthorList']
            full_author = None

        # Publication.select(graph, pm_id).first()
        publication = Publication.select(graph, pm_id).first()
        publication.title = title
        publication.volume = volume
        publication.issue = issue
        publication.pages = pages
        publication.year = date_of_pub
        publication.pubplace = pub_place
        publication.publisher = publisher
        graph.push(publication)
        create_author_nodes(publication, full_author)
    sys.stdout.write("Created Publications.")
    return pmid_set


def build_protein_interaction_rels(protein_interaction_dict):
    """
    Build protein-protein interactions
    :param protein_interaction_dict:
    :return:
    """
    for uni_id, interactors in protein_interaction_dict.items():
        if len(interactors) > 0:
            poly = Protein.select(graph, uni_id).first()
            interactors = interactors.split('; ')
            for interactor in interactors:
                if interactor == 'Itself':
                    interactor = poly.uniquename
                _poly = Protein.select(graph, interactor).first()
                if _poly:
                    poly.interacts_with.add(_poly)
                    graph.push(poly)


def create_drugbank_nodes():
    """
    Create DrugBank Drug Nodes
    :return:
    """
    sys.stdout.write("\nCreating DrugBank Nodes...\n")
    drug_set = set()
    zipped_tpi = zipfile.ZipFile(TARGET_PROTEIN_IDS)
    df = read_csv(zipped_tpi.open("all.csv")).fillna("")
    for entry in tqdm(df.values):
        uniprot_entry = entry[6]
        dbank_ids = entry[12]
        protein_ = Protein.select(graph).where(
            "_.entry_name='{}'".format(uniprot_entry))
        if protein_:
            for protein in protein_:
                drug_ids = [x.strip() for x in dbank_ids.split(';') if x]
                for _id in drug_ids:
                    drug_set.add(_id)
                    dbxref = DbXref(db="DrugBank", accession=_id)
                    graph.create(dbxref)
                    drug = Drug.select(graph, _id).first()
                    if drug:
                        drug.target.add(protein)
                    else:
                        drug = Drug(accession=_id)
                        graph.create(drug)
                    drug.target.add(protein)
                    graph.push(drug)
                    protein.drug.add(drug)
                    protein.dbxref.add(dbxref)
                    graph.push(protein)

    zipped_dv = zipfile.ZipFile(DRUG_VOCAB)
    df = read_csv(zipped_dv.open("drugbank vocabulary.csv")).fillna("")
    for entry in tqdm(df.values):
        dbank_id = entry[0].strip()
        commom_name = entry[2]
        synonyms = entry[5]
        drug = Drug.select(graph, dbank_id).first()
        if drug:
            drug.name = commom_name
            drug.synonyms = synonyms
            graph.push(drug)
    sys.stdout.write("\nDone Creating DrugBank Nodes...\n")


def map_gene_and_cds_to_protein(protein):
    """
    Mapping Proteins to CDS
    :param protein:
    :return:
    """
    if protein and protein.parent:
        for tag in protein.parent:
            gene = Gene.select(graph, tag).first()
            if not gene:
                gene = PseudoGene.select(graph, tag).first()
            if gene:
                gene.encodes.add(protein)
                graph.push(gene)
                protein.encoded_by.add(gene)
                graph.push(protein)
    # ens_id = map_ue_to_ens_trs(entry['Entry'])[0]
    protein_entry = str(protein.uniquename).strip()
    ens_id = eu_mapping(protein_entry, 'ENSEMBLGENOME_TRS_ID')
    if ens_id is not None:
        cds = CDS.select(graph, ens_id[0]).first()
        if cds:
            # Protein-[derives_from]->CDS
            protein.derives_from.add(cds)
            graph.push(protein)
            cds.derived.add(protein)
            graph.push(cds)


def split_gene_names(parent):
    """
    Split gene names that have spaces and /
    :param parent:
    :return:
    """
    if not parent:
        return None
    if ';' in parent:
        parent = parent.split(';')
    elif '/' in parent:
        parent = parent.split('/')
    else:
        parent = parent.split(' ')
    return parent


def create_protein_nodes():
    """
    Create Protein Nodes from UniProt results.
    :return:
    """
    sys.stdout.write("\nCreating Protein Nodes...\n")
    protein_interaction_dict = dict()

    df = read_csv(UNIPROT_DATA).fillna("")
    for entry in tqdm(df.values):
        protein_interaction_dict[str(entry[0]).strip()] = entry[6]
        dbxref = DbXref(db="UniProt", accession=entry[1], version=entry[0])
        graph.create(dbxref)

        # pdb_id = None
        # if entry[12] is not "":
        pdb_id = eu_mapping(entry[0], to='PDB_ID')
        protein = Protein()
        protein.name = entry[9]
        protein.uniquename = entry[0]
        protein.entry_name = entry[1]
        protein.ontology_id = protein.so_id
        protein.seqlen = entry[16]
        protein.residues = entry[14]
        # Catering for 'MT1076 MT1237 MT3197', 'MT0511/MT0512'
        # and 'Rv0132c LH57_00740'
        parent = split_gene_names(parent=entry[2])
        protein.parent = parent
        protein.family = entry[17]
        protein.function = entry[13]
        protein.pdb_id = pdb_id
        protein.mass = entry[15]
        protein.three_d = entry[12]
        graph.create(protein)
        protein.dbxref.add(dbxref)
        graph.push(protein)

        # create_chembl_nodes(protein, entry['Entry'])
        map_gene_and_cds_to_protein(protein)
        create_interpro_term_nodes(protein, entry[5])

    build_protein_interaction_rels(protein_interaction_dict)
    sys.stdout.write("\nCreated UniProt Nodes.")


def map_gene_to_orthologs(locus_tags):
    """
    Mapping Genes to orthologs
    :param locus_tags:
    :return:
    """
    sys.stdout.write("\nMapping Orthologs...\n")
    for tag_list in locus_tags:
        for tag in tag_list:
            gene = Gene.select(graph, tag).first()
            if gene:
                if tag.startswith('Rv'):
                    ortholog = fetch_ortholog(locus_tag=str(tag))
                    if ortholog:
                        orthologous_gene = Gene.select(
                            graph, str(ortholog)).first()
                        if orthologous_gene:
                            gene.orthologous_to.add(orthologous_gene)
                            orthologous_gene.orthologous_to_.add(gene)
                            graph.push(gene)
                            graph.push(orthologous_gene)
    sys.stdout.write("\nMapped Orthologs")


def create_kegg_pathways_nodes():
    """
    Create KEGG pathways
    :return:
    """
    sys.stdout.write("Creating KEGG Pathways...")

    def map_pathway_to_proteins(pathway_genes, path):
        for g_id in pathway_genes:
            g_id = "Rv" + \
                g_id.strip("RVBD_") if "RV" in g_id else g_id
            # Protein parent is stored as an array
            gene = Gene.select(graph, g_id).first()
            if gene:
                for protein in gene.encodes:
                    protein.pathway.add(path)
                    graph.push(protein)
                    path.protein.add(protein)
                    graph.push(path)

    # TODO: Add mtc
    organisms = ['mtu']
    for organism in organisms:
        kegg.organism = organism
        pathway_ids = kegg.pathwayIds
        for path in tqdm(pathway_ids):
            data = kegg.parse(kegg.get(path))
            accession = path[path.find(organism):].strip()
            if isinstance(data, dict) is True:
                pathway = Pathway()
                pathway.accession = accession
                pathway._class = data.get('CLASS')
                pathway.name = data['NAME'][0].replace(
                    " - Mycobacterium tuberculosis H37Rv", "")
                # pathway.name = data['PATHWAY_MAP'].get(path.strip("path:"))
                pathway.summation = data.get('DESCRIPTION')
                pathway.species = data.get('ORGANISM')
                graph.create(pathway)
                if data.get('GENE'):
                    map_pathway_to_proteins(data['GENE'].keys(), pathway)
            else:
                res_split = data.split("\n")

                # accession = res_split[0].split()[1]
                name = res_split[1].strip('NAME').strip().replace(
                    " - Mycobacterium tuberculosis H37Rv", "")
                summation = res_split[2].strip("DESCRIPTION").strip()
                _class = res_split[3].strip("CLASS").strip()
                # Create Pathway
                pathway = Pathway()
                pathway.accession = accession
                pathway.name = name
                pathway.summation = summation
                pathway._class = _class
                graph.create(pathway)

                gene_str_list = [
                    s.split() for s in res_split if
                    "Rv" in s and 'Myco' not in s
                ]
                genes = {
                    g for l in gene_str_list for g in l if
                    str(g).isalnum() and 'Rv' in g
                }

                map_pathway_to_proteins(genes, pathway)
    sys.stdout.write("\nCreated KEGG Pathway Nodes.")


def create_reactome_pathway_nodes():
    """
    Create REACTOME Pathway Nodes
    :return:
    """
    sys.stdout.write("\nCreating REACTOME Pathways...\n")
    with open(UNIPROT_DATA, 'r') as csv_file:
        reader = csv.DictReader(csv_file, delimiter=',')
        for entry in tqdm(reader):
            protein = entry['Entry']
            pathways = eu_mapping(protein, to='REACTOME_ID')
            if pathways:
                for pathway_id in pathways:
                    path_res = reactome_old.query_by_id("Pathway", pathway_id)
                    if not isinstance(path_res, int):
                        pathway = Pathway()
                        pathway.accession = pathway_id
                        pathway._type = path_res.get('schemaClass')
                        pathway.species = path_res.get('speciesName')
                        pathway.name = path_res.get('displayName', pathway_id)
                        if path_res.get('compartment'):
                            compartment = path_res.get('compartment')
                            pathway.compartment = compartment[0].get(
                                'displayName') if compartment else None
                        summation = path_res.get('summation')
                        pathway.summation = summation[0].get(
                            'displayName') if summation else None
                        graph.create(pathway)
                        _protein = Protein.select(graph, protein).first()
                        if _protein:
                            _protein.pathway.add(pathway)
                            graph.push(_protein)
                            pathway.protein.add(_protein)
                            graph.push(pathway)
    sys.stdout.write("\nCreated REACTOME Pathway Nodes.")


def create_known_mutation_nodes(**kwargs):
    """
    Create Known mutations
    :return:
    """
    fluoroquinolones = ["ciprofloxacin",
                        "ofloxacin", "levofloxacin", "moxifloxacin"]
    aminoglyconsides = ["amikacin", "kanamycin", "streptomycin", "capreomycin"]

    v_set = VariantSet(name=kwargs.get("vset_name", ""),
                       owner=kwargs.get("vset_owner", ""))
    call_set = CallSet(name=kwargs.get("cset_name", ""))

    v_set.has_callsets.add(call_set)
    call_set.belongs_to_vset.add(v_set)

    graph.create(v_set)
    graph.create(call_set)

    variant = Variant(chrom=kwargs.get("chrom", ""), pos=kwargs.get("pos", ""),
                      ref_allele=kwargs.get("ref_allele", ""),
                      alt_allele=kwargs.get("alt_allele", ""),
                      gene=kwargs.get("gene", ""), pk=kwargs.get("pk", ""),
                      consequence=kwargs.get("consequence", ""))
    variant.loc_in_seq = kwargs.get("loc_in_seq")
    variant.promoter = kwargs.get("promoter")
    variant.biotype = kwargs.get("biotype")
    variant.drug = kwargs.get("drug_name")
    variant.sources = kwargs.get("sources")
    variant.belongs_to_cset.add(call_set)
    call_set.has_variants.add(variant)

    def map_drug_class_to_variant(_class):
        """
        Map all drugs in class to variant
        :param _class:
        :return:
        """
        for item in _class:
            drugs = Drug.select(graph).where(
                "_.name=~'(?i).*{}.*'".format(item))
            for _drug in drugs:
                variant.resistant_to.add(_drug)

    if kwargs.get("drug_name") == "aminoglycosides":
        map_drug_class_to_variant(aminoglyconsides)
    elif kwargs.get("drug_name") == "fluoroquinolones":
        map_drug_class_to_variant(fluoroquinolones)
    else:
        for drug_id in kwargs.get("drugbank_id"):
            drug = Drug.select(graph, str(drug_id).upper()).first()
            if drug:
                variant.resistant_to.add(drug)
            elif drug_id and kwargs.get("drug_name"):
                drug = Drug(accession=drug_id.strip(),
                            name=kwargs.get("drug_name").capitalize())
                graph.create(drug)
                variant.resistant_to.add(drug)

    gene_name = str(kwargs.get("gene")).lower()
    gene = Gene.select(graph).where(
        f"_.name=~'(?i).*{gene_name}.*' OR _.uniquename=~'(?i).*{gene_name}.*'").first()
    if gene:
        variant.occurs_in.add(gene)
    else:
        rna = RRna.select(graph).where(
            f"_.name=~'(?i).*{gene_name}.*' OR _.uniquename=~'(?i).*{gene_name}.*'").first()
        if rna:
            variant.occurs_in_.add(rna)

    graph.create(variant)
    graph.push(call_set)


def create_operon_nodes(text_file=None):
    """
    Adding functional categories to Feature Nodes
    :param text_file:
    :return:
    """
    sys.stdout.write("\nAdding operon data...")
    with open(text_file) as text_file:
        for line in text_file:
            if 'OPERON' in str(line):
                tab_split = line.split('\t')
                # locus = tab_split[0]
                # gene_name = tab_split[1]
                # name_operon = tab_split[10]
                locus_operon = tab_split[11]
                description = tab_split[7]
                operon = Operon()
                # Must we use the product as the uniquename
                operon.uniquename = locus_operon
                operon.description = description
                graph.create(operon)
                genes = locus_operon.split(',')
                for locus_tag in genes:
                    gene = Gene.select(graph, locus_tag.strip()).first()
                    if gene:
                        gene.member_of.add(operon)
                        if len(genes) == 1:
                            gene.co_regulated.add(gene)
                        else:
                            # Let's not build reverse co-regulated rel
                            for g_id in genes[1:]:
                                g = Gene.select(graph, g_id.strip()).first()
                                if g:
                                    gene.co_regulated.add(g)
                        graph.push(gene)
                        operon.gene.add(gene)
                        graph.push(operon)


def map_srna_to_mrna(text_file):
    """
    Map sRNA to the mRNA they regulate
    :param text_file:
    :return:
    """
    sys.stdout.write("\nAdding sRNA data...")
    with open(text_file) as text_file:
        next(text_file)
        for line in text_file:
            tab_split = line.split('\t')
            srna_name = tab_split[0]
            # srna_fmax = tab_split[2]
            mrna_name = tab_split[7]
            # mrna_fmax = tab_split[9]
            ncrna = NCRna.select(graph).where(
                "_.name=~'(?i).*{}.*'".format(srna_name.lower())).first()
            if ncrna:
                mrnas = mrna_name.split()
                if "-" in mrna_name:
                    mrnas = mrna_name.split("-")
                for name in mrnas:
                    gene = Gene.select(graph).where(
                        "_.uniquename=~'(?i).*{}.*'".format(
                            name.lower())).first()
                    if gene:
                        ncrna.regulates_gene.add(gene)
                        graph.push(ncrna)
