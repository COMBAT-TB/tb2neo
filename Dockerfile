FROM python:3-slim
LABEL MAINTAINER="thoba@sanbi.ac.za"

RUN apt-get update -y --fix-missing && \
    apt-get upgrade -y && \
    apt-get install git -y && \
    mkdir /code && \
    pip install -U pip

COPY requirements.txt /code

RUN pip install -r /code/requirements.txt

COPY . /code
WORKDIR /code

RUN pip install -e .

ENTRYPOINT ["./docker-entrypoint.sh"]
