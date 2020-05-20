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
        git-lfs \
   && rm -rf /var/lib/apt/lists/*

WORKDIR /GNE

COPY . /GNE

RUN git clone https://github.com/kckishan/GNE.git

RUN pip install --trusted-host pypi.python.org -r requirements.txt

RUN ls ./

RUN ["python", "run_GNE.py"]







