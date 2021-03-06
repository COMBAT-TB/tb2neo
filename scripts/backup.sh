#!/usr/bin/env bash
# Back up named volumes created when deploying using docker-compose

NAMED_VOLUMES='tb2neo_db_data tb2neo_es_data tb2neo_loader_data'

for  VOL in ${NAMED_VOLUMES}
do
    VOLUME="$(docker volume ls -f name=${VOL} -q)"
    if [ ! -z ${VOLUME} ]; then
        echo "Backing up ${VOLUME} to ${HOME}..."
        docker run -v ${VOLUME}:/volume -v ${HOME}:/backup loomchild/volume-backup backup ${VOLUME}
    else
        echo "Unable to backup ${VOL}. It does not exist!"
    fi
done