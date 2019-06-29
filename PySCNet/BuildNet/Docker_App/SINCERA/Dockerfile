FROM r-base
LABEL maintainer = 'ming.wu@tum.de'

RUN apt-get update \
    && apt-get install -y \
        build-essential \
        libcurl4-openssl-dev \
        libssl-dev \
        curl \
        unzip \
        vim \
&& rm -rf /var/lib/apt/lists/*


WORKDIR /SINCERA

COPY . /SINCERA

RUN ["Rscript", "install_packages.r"]
