FROM mambaorg/micromamba:alpine 

LABEL org.opencontainers.image.source=https://github.com/metagenlab/assembly_finder
LABEL org.opencontainers.image.description="Snakemake-powered cli to download genomes with NCBI datasets"
LABEL org.opencontainers.image.licenses=MIT

RUN micromamba config prepend channels conda-forge && \
    micromamba config append channels bioconda && \
    micromamba config set channel_priority strict \
    micromamba config set extract_threads 1 && \
    micromamba install assembly_finder --only-deps -n base && \
    micromamba install mamba -n base && \
    micromamba clean -afy

COPY . .
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN python -m pip install -vvv --no-cache-dir --no-deps .
RUN assembly_finder -i bacteria -nb 1 --conda-create-envs-only --conda-cleanup-pkgs