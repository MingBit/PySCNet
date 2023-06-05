FROM julia:1.9.0

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

WORKDIR /PIDC

COPY . /PIDC

RUN ["julia", "run_pidc.jl"]



