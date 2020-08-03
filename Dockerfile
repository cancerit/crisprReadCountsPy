FROM ubuntu:20.04 as builder
USER root

RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
locales \
g++ \
make \
pkg-config \
python3 python3-dev python3-pip \
libbz2-dev liblzma-dev libcurl4-gnutls-dev
# libbz2-dev liblzma-dev libcurl4-gnutls-dev are for building pysam

ENV CGP_OPT /opt/wtsi-cgp
RUN mkdir $CGP_OPT
ENV PYTHONPATH $CGP_OPT/python-lib/lib/python3.8/site-packages

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

# build the tools in this repo, separate to reduce build time on errors
COPY . .
RUN python3 setup.py sdist
RUN pip3 install --prefix="$CGP_OPT/python-lib" dist/$(ls -1 dist/)

FROM ubuntu:20.04

LABEL maintainer="cgphelp@sanger.ac.uk" \
      uk.ac.sanger.cgp="Cancer, Ageing and Somatic Mutation, Wellcome Trust Sanger Institute" \
      version="2.0.0" \
      description="crisprReadCounts docker container"

RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
apt-transport-https \
locales \
ca-certificates \
time \
unattended-upgrades \
python3 \
libbz2-dev liblzma-dev libcurl4-gnutls-dev && \
unattended-upgrade -d -v && \
apt-get remove -yq unattended-upgrades && \
apt-get autoremove -yq

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

ENV CGP_OPT /opt/wtsi-cgp
ENV PATH $CGP_OPT/bin:$CGP_OPT/python-lib/bin:$PATH
ENV PYTHONPATH $CGP_OPT/python-lib/lib/python3.8/site-packages
ENV LD_LIBRARY_PATH $OPT/lib
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

RUN mkdir -p $CGP_OPT
COPY --from=builder $CGP_OPT $CGP_OPT

## USER CONFIGURATION
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER ubuntu
WORKDIR /home/ubuntu
# Add a VERSION file into the container, jsut another way to veriry which version of code is in this container
COPY crispr_read_counts/version.py .
RUN echo $(grep 'version=' version.py | cut -d "'" -f 2) > VERSION && rm version.py

CMD ["/bin/bash"]
