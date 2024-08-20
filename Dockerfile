FROM mambaorg/micromamba:1.5.8
LABEL org.opencontainers.image.source=https://github.com/metagenlab/assembly_finder
LABEL org.opencontainers.image.description="Snakemake-powered cli to download genomes with NCBI datasets"
LABEL org.opencontainers.image.licenses=MIT
ENV LANG=C.UTF-8 
ENV SHELL=/bin/bash 
USER root

COPY . /pkg
RUN micromamba config set extract_threads 1 && \
    micromamba create --always-copy -p /env -y -f /pkg/env.lock && \
    micromamba clean -afy 

RUN micromamba run -p /env python -m pip install /pkg --no-deps --no-build-isolation --no-cache-dir -vvv
ENV PATH="/env/bin:$PATH" XDG_CACHE_HOME=/tmp