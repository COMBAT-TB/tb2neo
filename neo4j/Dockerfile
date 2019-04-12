FROM neo4j:3.4
LABEL Maintainer="lose.thoba@gmail.com"

ENV GRAPHAWARE_PROD_URL https://products.graphaware.com/download

RUN apk update \
    && apk add --no-cache ca-certificates wget curl
#     && wget ${GRAPHAWARE_PROD_URL}/framework-server-community/graphaware-server-community-all-3.3.0.51.jar \
#     -P plugins/ \
#     && wget ${GRAPHAWARE_PROD_URL}/neo4j-to-elasticsearch/graphaware-neo4j-to-elasticsearch-3.3.0.51.7.jar \
#     -P plugins/ \
#     && wget ${GRAPHAWARE_PROD_URL}/uuid/graphaware-uuid-3.3.0.51.14.jar \
#     -P plugins/ \
#     && ls plugins/
RUN wget https://github.com/neo4j-contrib/neo4j-elasticsearch/releases/download/3.4.5/neo4j-elasticsearch-3.4.5.jar \
    -P plugins/ \
    && wget https://github.com/neo4j-contrib/neo4j-graph-algorithms/releases/download/3.4.12.7/graph-algorithms-algo-3.4.12.7.jar \
    -P plugins/ \
    && wget https://github.com/neo4j-contrib/neo4j-apoc-procedures/releases/download/3.4.0.5/apoc-3.4.0.5-all.jar \
    -P plugins/ 
COPY neo4j.conf conf/

ENV  NEO4J_dbms_security_procedures_unrestricted='apoc.*,algo.*'
RUN echo 'browser.remote_content_hostname_whitelist=*' >> conf/neo4j.conf


