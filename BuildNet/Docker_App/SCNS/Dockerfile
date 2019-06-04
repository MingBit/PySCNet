FROM fsharp

LABEL maintainer='ming.wu@tum.de'


RUN apt-get update \
    && apt-get install -y \
        git \
        curl \
        unzip \
        vim \
        wget \
   && rm -rf /var/lib/apt/lists/*

WORKDIR /SCNS

COPY . /SCNS

RUN  git clone https://github.com/swoodhouse/SCNS-Toolkit.git

RUN ls ./




