FROM python:2.7-slim
MAINTAINER Thoba Lose "thoba@sanbi.ac.za"

RUN apt-get update -y --fix-missing && apt-get upgrade -y
RUN mkdir /code && \
    pip install -U pip
COPY requirements.txt /code
RUN apt-get install git -y && \
    pip install git+https://github.com/cokelaer/bioservices.git
RUN pip install -r /code/requirements.txt
COPY . /code
WORKDIR /code
CMD ["python", "main.py"]