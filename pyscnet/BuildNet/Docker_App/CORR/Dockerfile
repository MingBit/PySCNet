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


WORKDIR /CORR
COPY . /CORR

RUN Rscript -e "install.packages('ppcor')"
RUN Rscript -e "install.packages('lsa')"
RUN ["Rscript", "run_corr.r", "/CORR/Expr.txt"]



