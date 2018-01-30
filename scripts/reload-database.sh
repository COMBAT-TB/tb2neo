#!/usr/bin/env bash
# Checks if the gff2neo container exists and has exited. It then restarts it if it does.
# Else run a container from the gff2neo image

CONT_NAME="gff2neo_gff2neo_1"
IMAGE_NAME="gff2neo_gff2neo"
CID="$(docker ps -aq -f status=exited -f name=${CONT_NAME})"

if [ CID ]; then
 echo "Restarting ${CONT_NAME} with ID: ${CID}..."
 docker restart ${CONT_NAME}
else
 echo "Running ${CONT_NAME} from ${IMAGE_NAME}..."
 docker run -d --name ${CONT_NAME} ${IMAGE_NAME}
fi
