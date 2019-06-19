FROM neo4j:3.5-enterprise
LABEL Maintainer="lose.thoba@gmail.com"

RUN apt-get update \
    && apt-get install -y wget openssl apt-utils

RUN wget https://github.com/neo4j-contrib/neo4j-elasticsearch/releases/download/3.5.6/neo4j-elasticsearch-3.5.6.jar \
    -P plugins/ \
    && wget https://github.com/neo4j-contrib/neo4j-graph-algorithms/releases/download/3.5.4.0/graph-algorithms-algo-3.5.4.0.jar \
    -P plugins/ \
    && wget https://github.com/neo4j-contrib/neo4j-apoc-procedures/releases/download/3.4.0.7/apoc-3.4.0.7-all.jar \
    -P plugins/ 

# COPY neo4j.conf conf/

ENV NEO4J_AUTH=none \
    NEO4J_ACCEPT_LICENSE_AGREEMENT=yes

# APOC and Algorithms settings
RUN echo 'dbms.security.procedures.whitelist=apoc.*, algo.*' >> conf/neo4j.conf \
    && echo 'dbms.security.procedures.unrestricted=apoc.*, algo.*' >> conf/neo4j.conf \
    && echo 'browser.remote_content_hostname_whitelist=*' >> conf/neo4j.conf \
    # Neo4j Elastic{Search} Integration
    && echo 'elasticsearch.host_name=http://localhost:9200'  >> conf/neo4j.conf \
    && echo 'elasticsearch.index_spec=gene:Gene(uniquename, name, description, biotype, category), mrna:MRna(uniquename, name, parent, biotype, category), transcript:Transcript(uniquename, name, parent, biotype, category), trna:TRna(uniquename, name, parent, biotype, category), ncrna:NCRna(uniquename, name, parent, biotype, category), rrna:RRna(uniquename, name, parent, biotype, category), pseudogene:PseudoGene(uniquename, name, description, biotype, category), protein:Protein(uniquename, name, function, family, parent), pathway:Pathway(accession, name, summation), drug:Drug(accession, name, synonyms, definition), goterm:GOTerm(name, accession, definition), publication:Publication(title), location:Location(fmin, fmax), variant:Variant(pos, impact, gene, gene_id, known, consequence, drug)' >> conf/neo4j.conf \
    && echo 'elasticsearch.include_id_field=false' >> conf/neo4j.conf \
    && echo 'elasticsearch.include_labels_field=false' >> conf/neo4j.conf \
    && echo 'elasticsearch.discovery=true' >> conf/neo4j.conf



