"""
Interface to the Neo4j Database
"""
import os
import sys

from bioservices import ChEMBL, QuickGO, Reactome, KEGG
from py2neo import Graph

from model.core import *
from ncbi import fetch_publication_list
from quickgo import fetch_quick_go_data
from uniprot import *

graph = Graph(host=os.environ.get("DB", "localhost"), bolt=True,
              password=os.environ.get("NEO4J_PASSWORD", ""))

chembl = ChEMBL(verbose=False)
quick_go = QuickGO(verbose=False)
reactome = Reactome(verbose=False)
kegg = KEGG(verbose=False)

# watch("neo4j.bolt")

gene_dict = dict()
transcript_dict = dict()
pseudogene_dict = dict()
cds_dict = dict()
exon_dict = dict()
rrna_dict = dict()
trna_dict = dict()
ncrna_dict = dict()
location_dict = dict()
go_term_set = set()

target_protein_ids_csv = "data/drugbank/all_target_polypeptide_ids.csv"
drug_vocab_csv = "data/drugbank/drugbank_vocabulary.csv"


def delete_data():
    """
    Delete existing data.
    :return:
    """
    # print("Deleting all nodes and relationships in {}".format(graph))
    sys.stdout.write(
        "Deleting all nodes and relationships in {}".format(graph))

    graph.delete_all()


def create_organism_nodes():
    """
    Create Organism Nodes
    :return:
    """
    abbrev = "H37Rv"
    strain = abbrev
    genus = "Mycobacterium"
    species = "M. tuberculosis"
    common_name = "TB"

    organism = Organism(abbreviation=abbrev, strain=strain, genus=genus,
                        species=species, common_name=common_name)
    graph.create(organism)


def create_chromosome_nodes():
    """
    Create Chromosome Nodes
    :return:
    """
    name = "Chr1"
    uniquename = "Chr1"
    chromosome = Chromosome()
    chromosome.name = name
    chromosome.uniquename = uniquename
    graph.create(chromosome)


def create_gene_nodes(feature):
    """
    Create Gene Nodes
    :param feature:
    :return:
    """
    names = get_feature_name(feature)
    name = names.get("Name", names.get("UniqueName"))
    unique_name = names.get("UniqueName", name)
    description = feature.qualifiers["description"]
    biotype = feature.qualifiers['biotype'][0]
    parent = get_feature_parent(feature)

    gene = Gene()
    gene.name = name
    gene.uniquename = unique_name
    gene.parent = parent
    gene.biotype = biotype
    gene.description = description
    graph.create(gene)
    gene_dict[unique_name] = gene


def create_transcript_nodes(feature):
    """
    Create Transcipt Nodes
    :param feature:
    :return:
    """
    names = get_feature_name(feature)
    name = names.get("Name", names.get("UniqueName"))
    unique_name = names.get("UniqueName", name)
    biotype = feature.qualifiers['biotype'][0]
    parent = get_feature_parent(feature)

    transcript = Transcript()
    transcript.name = name
    transcript.parent = parent
    transcript.uniquename = unique_name
    transcript.biotype = biotype
    graph.create(transcript)
    transcript_dict[unique_name] = transcript


def create_pseudogene_nodes(feature):
    """
    Create Pseudogene Nodes
    :param feature:
    :return:
    """
    names = get_feature_name(feature)
    name = names.get("Name", names.get("UniqueName"))
    unique_name = names.get("UniqueName", name)
    description = feature.qualifiers["description"][0]
    biotype = feature.qualifiers['biotype'][0]
    parent = get_feature_parent(feature)

    pseudogene = PseudoGene()
    pseudogene.name = name
    pseudogene.uniquename = unique_name
    pseudogene.parent = parent
    pseudogene.description = description
    pseudogene.biotype = biotype
    graph.create(pseudogene)
    pseudogene_dict[unique_name] = pseudogene


def create_exon_nodes(feature):
    """
    Create Exon Nodes
    :param feature:
    :return:
    """
    names = get_feature_name(feature)
    name = names.get("Name", names.get("UniqueName"))
    unique_name = names.get("UniqueName", name)
    parent = get_feature_parent(feature)

    exon = Exon()
    exon.name = name
    exon.uniquename = unique_name
    exon.parent = parent
    graph.create(exon)
    exon_dict[unique_name] = exon


def create_rna_nodes(feature):
    """
    Create RNA Nodes
    :param feature:
    :return:
    """
    names = get_feature_name(feature)
    name = names.get("Name", names.get("UniqueName"))
    unique_name = names.get("UniqueName", name)
    parent = get_feature_parent(feature)

    if feature.type == 'tRNA_gene':
        trna = TRna()
        trna.name = name
        trna.parent = parent
        trna.uniquename = unique_name
        graph.create(trna)
        trna_dict[unique_name] = trna
    if feature.type == 'ncRNA_gene':
        ncrna = NCRna()
        ncrna.name = name
        ncrna.parent = parent
        ncrna.uniquename = unique_name
        graph.create(ncrna)
        ncrna_dict[unique_name] = ncrna
    if feature.type == 'rRNA_gene':
        rrna = RRna()
        rrna.name = name
        rrna.parent = parent
        rrna.uniquename = unique_name
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
    parent = get_feature_parent(feature)

    cds = CDS()
    cds.name = name
    cds.parent = parent
    cds.uniquename = unique_name
    graph.create(cds)
    cds_dict[unique_name] = cds


def get_feature_parent(feature):
    if feature.qualifiers.get('Parent'):
        parent = feature.qualifiers['Parent'][0]
        # [feature.qualifiers['Parent'][0].find(":") + 1:]
    else:
        parent = None
    return parent


def create_featureloc_nodes(feature):
    """
    Create FeatureLoc Nodes
    :param feature:
    :return:
    """
    srcfeature_id = get_feature_name(feature).get("UniqueName")
    pk = feature.location.start + feature.location.end
    feature_loc = Location(pk=pk, fmin=feature.location.start, fmax=feature.location.end,
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
    if feature.qualifiers.get("Name"):
        names["Name"] = feature.qualifiers["Name"][0]
        names["UniqueName"] = feature.id
        # [feature.id.find(":") + 1:]
    else:
        names["Name"] = names["UniqueName"] = feature.id
        # [feature.id.find(":") + 1:]
    return names


def build_gff_rels():
    """
    Build GFF Feature relationships
    :return:
    """
    sys.stdout.write("\nBuilding GFF Relationships...\n")
    for t, transcript in transcript_dict.iteritems():
        if transcript.parent in gene_dict.keys():
            gene = gene_dict.get(transcript.parent)
            transcript.part_of_g.add(gene)
            graph.push(transcript)
            gene.part_of.add(transcript)
            graph.push(gene)
        for p, pseudogene in pseudogene_dict.iteritems():
            if transcript.parent == pseudogene.uniquename:
                transcript.part_of_pg.add(pseudogene)
                graph.push(transcript)
                pseudogene.part_of.add(transcript)
                graph.push(pseudogene)
        for c, cds in cds_dict.iteritems():
            if transcript.uniquename == cds.parent:
                cds.part_of.add(transcript)
                graph.push(cds)
                transcript.part_of_cds.add(cds)
                graph.push(transcript)


def map_to_location(feature):
    """
    Map GFF Features to Location and Organism
    :param feature:
    :return:
    """
    # Find feature location with a srcfeature_id attr. matching this features uniquename and link them via
    # LOCATED_AT
    srcfeature_id = get_feature_name(feature).get("UniqueName")
    location = location_dict.get(srcfeature_id)
    organism = Organism.select(graph).first()
    chromosome = Chromosome.select(graph).first()
    rna = ["tRNA_gene", "ncRNA_gene", "rRNA_gene"]
    if location:
        if feature.type == 'gene':
            _feature = gene_dict.get(srcfeature_id)
            _feature.location.add(location)
            _feature.belongs_to.add(organism)
            _feature.located_on.add(chromosome)
            graph.push(_feature)
        elif feature.type == 'pseudogene':
            _feature = pseudogene_dict.get(srcfeature_id)
            _feature.location.add(location)
            _feature.belongs_to.add(organism)
            _feature.located_on.add(chromosome)
            graph.push(_feature)
        elif feature.type == 'exon':
            _feature = exon_dict.get(srcfeature_id)
            _feature.location.add(location)
            _feature.belongs_to.add(organism)
            _feature.located_on.add(chromosome)
            graph.push(_feature)
        elif feature.type in rna:
            if feature.type == 'tRNA_gene':
                _feature = trna_dict.get(srcfeature_id)
                _feature.location.add(location)
                _feature.belongs_to.add(organism)
                _feature.located_on.add(chromosome)
                graph.push(_feature)
            if feature.type == 'ncRNA_gene':
                _feature = ncrna_dict.get(srcfeature_id)
                _feature.location.add(location)
                _feature.belongs_to.add(organism)
                _feature.located_on.add(chromosome)
                graph.push(_feature)
            if feature.type == 'rRNA_gene':
                _feature = rrna_dict.get(srcfeature_id)
                _feature.location.add(location)
                _feature.belongs_to.add(organism)
                _feature.located_on.add(chromosome)
                graph.push(_feature)
        elif feature.type == 'CDS':
            _feature = cds_dict.get(srcfeature_id)
            _feature.location.add(location)
            _feature.belongs_to.add(organism)
            _feature.located_on.add(chromosome)
            graph.push(_feature)
        elif feature.type == 'transcript':
            _feature = transcript_dict.get(srcfeature_id)
            _feature.location.add(location)
            _feature.belongs_to.add(organism)
            _feature.located_on.add(chromosome)
            graph.push(_feature)


def create_is_a_cv_term_rel(go_set):
    """
    Creating IS_A relationships between CVTerms
    :param go_set: set of GO ids
    :return:
    """
    for go_id in go_set:
        is_a_list = fetch_quick_go_data(go_id)
        go_term = GOTerm.select(graph, go_id).first()
        for go in is_a_list:
            goid = go[go.find('G'):go.find('!')].strip()
            term = GOTerm.select(graph, goid).first()
            if term and go_term:
                go_term.is_a.add(term)
                graph.push(go_term)


def create_go_term_nodes():
    """
    Create CvTerm Nodes and build Polypetide relationships.
    :param protein:
    :param bp:
    :param cc:
    :param mf:
    :return:
    """
    sys.stdout.write("\nCreating GoTerm Nodes...\n")
    with open(uniprot_data_csv, 'rb') as csv_file:
        import time
        start = time.time()
        reader = csv.DictReader(csv_file, delimiter=',')
        for entry in reader:
            protein_entry = entry['Entry']
            protein = None
            if protein_entry is not '':
                protein = Protein.select(graph, protein_entry).first()
            go_ids = [g for g in entry['GO_IDs'].split("; ") if g is not '']
            for go_id in go_ids:
                # will map is_a
                go_term_set.add(go_id)
                result = quick_go.Term(go_id, frmt="obo").split('\n')
                name = result[2].split(":")[1]
                _def = result[3].split(":")[1]
                go_term = GOTerm(accession=go_id, name=name.strip(), definition=_def.strip())
                graph.create(go_term)
                if protein is not None:
                    protein.assoc_goterm.add(go_term)
                    graph.push(protein)
                    go_term.protein.add(protein)
                    graph.push(go_term)
        end = time.time()
        print("Created {} in {} seconds.".format(len(go_term_set), end - start))
    create_is_a_cv_term_rel(go_term_set)


def create_interpro_term_nodes(protein, entry):
    """
    Create InterPro Term Nodes.
    :param protein:
    :param entry:
    :return:
    """
    # http://generic-model-organism-system-database.450254.n5.nabble.com/Re-GMOD-devel-Storing-Interpro-domains-in-Chado-td459778.html
    terms = [t for t in entry.split("; ") if t is not '']
    for interpro in terms:
        import time
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


def create_pub_nodes():
    p_id_set = set()
    sys.stdout.write("\nCreating Publication Nodes...\n")
    import time
    _start = time.time()
    with open(uniprot_data_csv, 'rb') as csv_file:

        reader = csv.DictReader(csv_file, delimiter=',')
        for entry in reader:
            protein_entry = entry['Entry']
            protein = None
            if protein_entry is not '':
                protein = Protein.select(graph, protein_entry).first()
            pubmed_ids = [p for p in entry['PubMed'].split("; ") if p is not '']
            for p_id in set(pubmed_ids):
                pub = Publication()
                pub.pmid = p_id
                graph.create(pub)
                p_id_set.add(p_id)
                if protein:
                    protein.published_in.add(pub)
                    graph.push(protein)
    # Let's build a set before calling API and DB
    num_ids = len(p_id_set)
    chunksize = 500
    records = []
    for start in range(0, num_ids, chunksize):
        subset = list(p_id_set)[start:start + chunksize]
        records.extend(fetch_publication_list(subset))
    record_loaded_count = 0
    for record in records:
        if len(record) < 2:
            pm_id = record['id:'][0][record['id:'][0].find('able: ') + 6:]
            record = fetch_publication_list(pm_id, rettype='xml')
            rec = next(record)
            article = rec['MedlineCitation']['Article']
            title = article['ArticleTitle']
            pages = article['Pagination']['MedlinePgn']
            volume = article['Journal']['JournalIssue']['Volume']
            issue = article['Journal']['JournalIssue']['Issue']
            date_of_pub = article['Journal']['JournalIssue']['PubDate']['Month'] + " " + \
                          article['Journal']['JournalIssue']['PubDate']['Year']
            pub_place = rec['MedlineCitation']['MedlineJournalInfo']['Country']
            publisher = None
            author = None
            # full_author = article['AuthorList']
            full_author = None
        else:
            # https://www.nlm.nih.gov/bsd/mms/medlineelements.html
            pm_id = record['PMID']
            # there is internal caching so using a dictionary here doesn't
            # actually seem to save any time - pvh
            title = record.get('TI', None)
            volume = record.get('VI', None)
            issue = record.get('IP', None)
            pages = record.get('PG', None)
            date_of_pub = record.get('DP', None)
            pub_place = record.get('PL', None)
            publisher = record.get('SO', None)
            author = record.get('AU', None)
            full_author = record.get('FAU', None)

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
        record_loaded_count += 1
    end = time.time()
    sys.stdout.write("Created Publications in {} seconds.".format(end - _start))


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
                if _poly is None:
                    print("No Protein with uniquename: {}".format(interactor))
                    # time.sleep(2)
                else:
                    poly.interacts_with.add(_poly)
                    graph.push(poly)


def create_chembl_nodes(protein, entry):
    """
    Create ChEMBL Drug nodes from UniProt results.
    :return:
    """
    target = chembl.get_target_by_uniprotId(entry)
    if not isinstance(target, int):
        dbxref = DbXref(db="ChEMBL", accession=target['chemblId'])
        graph.create(dbxref)
        drug = Drug(accession=target['chemblId'], name=target["preferredName"], synonyms=target["synonyms"],
                    definition=target["description"])
        graph.create(drug)
        drug.target.add(protein)
        graph.push(drug)
        protein.drug.add(drug)
        protein.dbxref.add(dbxref)
        graph.push(protein)


def create_drugbank_nodes():
    """
    Create DrugBank Drug Nodes
    :return:
    """
    sys.stdout.write("\nCreating DrugBank Nodes...\n")
    drug_set = set()
    with open(target_protein_ids_csv, "rb") as csv_file:
        reader = csv.DictReader(csv_file)
        for _target in reader:
            if 'tuberculosis' in _target['Species']:
                print(_target['Gene Name'], _target['Uniprot Title'], _target['Drug IDs'])
                protein_ = Protein.select(graph).where("_.entry_name='{}'".format(_target['Uniprot Title']))
                if protein_:
                    for protein in protein_:
                        drug_ids = [x for x in _target['Drug IDs'].split('; ') if x is not '']
                        for _id in drug_ids:
                            drug_set.add(_id)
                            dbxref = DbXref(db="DrugBank", accession=_id)
                            graph.create(dbxref)
                            drug = Drug(accession=_id)
                            graph.create(drug)
                            drug.target.add(protein)
                            graph.push(drug)
                            protein.drug.add(drug)
                            protein.dbxref.add(dbxref)
                            graph.push(protein)

    with open(drug_vocab_csv, "rb") as csv_file_:
        reader = csv.DictReader(csv_file_)
        for entry in reader:
            for drug_id in drug_set:
                if entry['DrugBank ID'] == drug_id:
                    drug = Drug.select(graph, drug_id).first()
                    drug.name = entry['Common name']
                    drug.synonyms = entry['Synonyms']
                    graph.push(drug)


def map_cds_to_protein(protein, entry):
    """
    Mapping Proteins to CDS
    :param protein:
    :param entry:
    :return:
    """
    # Map CDS to Protein
    # ens_id = map_ue_to_ens_trs(entry['Entry'])[0]
    ens_id = eu_mapping(entry, to='ENSEMBLGENOME_TRS_ID')
    if ens_id is not None and 'CCP' in ens_id[0]:
        cds = CDS.select(graph, "CDS:" + ens_id[0]).first()
        if cds:
            # Polypetide-derives_from->CDS
            protein.derives_from.add(cds)
            graph.push(protein)
            cds.derived.add(protein)
            graph.push(cds)


def create_uniprot_nodes():
    """
    Build DbXref nodes from UniProt results.
    :return:
    """
    sys.stdout.write("\nCreating UniProt Nodes from CSV...\n")
    start = time()
    protein_interaction_dict = dict()
    with open(uniprot_data_csv, 'rb') as csv_file:
        reader = csv.DictReader(csv_file, delimiter=',')
        for entry in reader:
            protein_interaction_dict[entry['Entry']] = entry['Interacts_With']
            dbxref = DbXref(db="UniProt", accession=entry[
                'Entry_Name'], version=entry['Entry'])
            graph.create(dbxref)

            pdb_id = None
            if len(entry['3D']) > 0:
                pdb_id = eu_mapping(entry['Entry'], to='PDB_ID')
            protein = Protein()
            protein.name = entry['Protein_Names']
            protein.uniquename = entry['Entry']
            protein.entry_name = entry['Entry_Name']
            protein.ontology_id = protein.so_id
            protein.seqlen = entry['Length']
            protein.residues = entry['Sequence']
            protein.parent = entry['Gene_Names_OL']
            protein.family = entry['Protein_Families']
            protein.function = entry['Function_CC']
            protein.pdb_id = pdb_id
            protein.mass = entry['Mass']
            protein.three_d = entry['3D']
            graph.create(protein)
            protein.dbxref.add(dbxref)
            graph.push(protein)

            # create_chembl_nodes(protein, entry['Entry'])
            map_cds_to_protein(protein, entry['Entry'])

            create_interpro_term_nodes(protein, entry['InterPro'])
    build_protein_interaction_rels(protein_interaction_dict)
    end = time()
    print("\nDone creating UniProt Nodes in ", end - start, "secs.")


def map_gene_protein(locus_tags):
    sys.stdout.write("\nMapping Genes to Proteins...\n")
    start = time()
    for tag_list in locus_tags:
        tags = ["gene:" + x for x in tag_list]
        for tag in tags:
            gene = Gene.select(graph, tag).first()
            if gene:
                parent = gene.uniquename[gene.uniquename.find(':') + 1:]
                protein = Protein.select(graph).where("_.parent='{}'".format(parent)).first()
                if protein:
                    gene.encodes.add(protein)
                    graph.push(gene)
                    protein.encoded_by.add(gene)
                    graph.push(protein)
                print(gene.uniquename, [protein.uniquename if protein else None])
    end = time()
    print("\nDone mapping Genes to Proteins in ", end - start, "secs.")


def create_kegg_pathways():
    sys.stdout.write("Creating KEGG Pathways...")
    start = time()
    kegg.organism = 'mtv'
    pathway_ids = kegg.pathwayIds
    for path in pathway_ids:
        data = kegg.parse(kegg.get(path))
        pathway = Pathway()
        pathway.accession = path[path.find('mtv'):].strip()
        pathway._class = data.get('CLASS')
        pathway.name = data['PATHWAY_MAP'].get(path.strip("path:"))
        pathway.summation = data.get('DESCRIPTION')
        pathway.species = data.get('ORGANISM')
        graph.create(pathway)
        if data.get('GENE'):
            for g_id in data['GENE'].keys():
                protein_ = Protein.select(graph).where("_.parent='{}'".format("Rv" + g_id.strip("RVBD_")))
                if protein_:
                    for protein in protein_:
                        protein.pathway.add(pathway)
                        graph.push(protein)
                        pathway.protein.add(protein)
                        graph.push(pathway)
    end = time()
    sys.stdout.write("\nDone creating KEGG Pathway Nodes in {} secs.".format(end - start))


def create_pathway_nodes():
    """
    Create Pathway Nodes
    :return:
    """
    sys.stdout.write("\nCreating REACTOME Pathways...\n")
    start = time()
    with open(uniprot_data_csv, 'rb') as csv_file:
        reader = csv.DictReader(csv_file, delimiter=',')
        for entry in reader:
            protein = entry['Entry']
            pathways = eu_mapping(protein, to='REACTOME_ID')
            if pathways:
                for pathway_id in pathways:
                    path_res = reactome.query_by_id("Pathway", pathway_id)
                    if not isinstance(path_res, int):
                        print(protein, "Pathway: {} - {}".format(pathway_id, path_res['displayName']))
                        pathway = Pathway()
                        pathway.accession = pathway_id
                        pathway._type = path_res['schemaClass']
                        pathway.species = path_res['speciesName']
                        pathway.name = path_res['displayName']
                        if path_res.get('compartment'):
                            pathway.compartment = path_res['compartment'][0]['displayName']
                        pathway.summation = path_res['summation'][0]['displayName']
                        graph.create(pathway)
                        _protein = Protein.select(graph, protein).first()
                        if _protein:
                            _protein.pathway.add(pathway)
                            graph.push(_protein)
                            pathway.protein.add(_protein)
                            graph.push(pathway)
    end = time()
    sys.stdout.write("\nDone creating REACTOME Pathway Nodes in {} secs.".format(end - start))
