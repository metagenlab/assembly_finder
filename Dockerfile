FROM mambaorg/micromamba
LABEL org.opencontainers.image.source=https://github.com/metagenlab/assembly_finder
LABEL org.opencontainers.image.description="Snakemake-powered cli to download genomes with NCBI datasets"
LABEL org.opencontainers.image.licenses=MIT
ENV LANG=C.UTF-8 SHELL=/bin/bash PATH="$MAMBA_ROOT_PREFIX/bin:$PATH" XDG_CACHE_HOME=/tmp

COPY --chown=$MAMBA_USER:$MAMBA_USER . /tmp
RUN micromamba config set extract_threads 1 && \
    micromamba install -n base -y -f /tmp/env.lock && \
    micromamba clean -afy 

ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN python -m pip install /tmp --no-deps --no-build-isolation --no-cache-dir -vvv