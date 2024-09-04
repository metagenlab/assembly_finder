FROM mambaorg/micromamba
LABEL org.opencontainers.image.source=https://github.com/metagenlab/assembly_finder
LABEL org.opencontainers.image.description="Snakemake-powered cli to download genomes with NCBI datasets"
LABEL org.opencontainers.image.licenses=MIT
ENV LANG=C.UTF-8 
ENV SHELL=/bin/bash

USER root
ENV APT_PKGS="curl rsync p7zip-full procps ca-certificates"
RUN apt-get update \
    && apt-get install -y --no-install-recommends ${APT_PKGS} \
    && apt-get clean \
    && rm -rf /var/lib/apt /var/lib/dpkg /var/lib/cache /var/lib/log
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER . /pkg
RUN micromamba config set extract_threads 1 && \
    micromamba install -n base -y -f /pkg/env.yaml && \
    micromamba clean -afy 
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install /pkg --no-deps --no-build-isolation --no-cache-dir -vvv
ENV PATH="/opt/conda/bin:$PATH" XDG_CACHE_HOME=/tmp
