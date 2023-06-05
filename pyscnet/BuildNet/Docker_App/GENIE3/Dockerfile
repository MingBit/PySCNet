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

WORKDIR /GENIE3

COPY . /GENIE3

RUN pip install --upgrade pip

RUN pip install --trusted-host pypi.python.org -r requirements

RUN ls ./

RUN ["python", "run_genie3.py"]






