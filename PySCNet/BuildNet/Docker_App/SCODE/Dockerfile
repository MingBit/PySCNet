FROM python:3.7

LABEL maintainer='ming.wu@tum.de'


RUN apt-get update \
    && apt-get install -y \
        ruby ruby-dev \
        r-base \
        build-essential \
        cmake \
        git \
        curl \
        unzip \
        vim \
   && rm -rf /var/lib/apt/lists/*



RUN git clone https://github.com/hmatsu1226/SCODE
COPY . /SCODE

RUN pip install pandas

WORKDIR SCODE/
RUN ["python", "run_scode.py"]

