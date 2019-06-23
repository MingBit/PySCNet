FROM python:3.6.1
LABEL maintainer = 'ming.wu@tum.de'

USER root
WORKDIR /SJARACNE
COPY . /SJARACNE

RUN apt-get install git
RUN git clone https://github.com/jyyulab/SJARACNe
RUN pip install -r requirements.txt
RUN mv setup.py SJARACNe/






