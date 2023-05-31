FROM python:3.7

LABEL maintainer='ming.wu@tum.de'


RUN apt-get update \
    && apt-get install -y \
        build-essential \
        cmake \
        git \
        curl \
        unzip \
        vim \
   && rm -rf /var/lib/apt/lists/*

WORKDIR /SCNODE2VEC

COPY . /SCNODE2VEC

RUN curl https://raw.githubusercontent.com/MingBit/SCNode2Vec/master/src/model_node2vec.py > scnode2vec_raw.py

RUN pip install --upgrade pip

RUN pip install --trusted-host pypi.python.org -r requirements.txt

RUN ls ./

RUN ["python", "run_scnode2vec.py"]







