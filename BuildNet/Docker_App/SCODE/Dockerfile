FROM r-base
LABEL maintainer = 'ming.wu@tum.de'

RUN apt-get update \
    && apt-get install -y \
        build-essential \
        libcurl4-openssl-dev \
        libssl-dev \
        curl \
        git \
        vim \
&& rm -rf /var/lib/apt/lists/*


WORKDIR /SCODE
COPY . /SCODE

RUN git clone https://github.com/hmatsu1226/SCODE
RUN cd SCODE

