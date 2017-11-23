FROM neo4j:3.2.7
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
RUN wget https://github.com/neo4j-contrib/neo4j-elasticsearch/releases/download/3.2.3/neo4j-elasticsearch-3.2.3.jar \
    -P plugins/
COPY neo4j.conf conf/
# create es index
# RUN curl curl -XPUT 'http://es:9200/neo4j-index' \
#     && apk del wget ca-certificates


