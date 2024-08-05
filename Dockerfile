FROM mambaorg/micromamba
LABEL org.opencontainers.image.source=https://github.com/metagenlab/assembly_finder
LABEL org.opencontainers.image.description="Snakemake-powered cli to download genomes with NCBI datasets"
LABEL org.opencontainers.image.licenses=MIT
ADD . /tmp/repo
WORKDIR /tmp/repo
ENV LANG C.UTF-8
ENV SHELL /bin/bash
USER root 

RUN micromamba config set extract_threads 1 && \
    micromamba install -q -y -c bioconda -c conda-forge -n base \
    assembly_finder --only-deps && \
    micromamba install -q -y -c conda-forge -n base mamba && \
    micromamba clean -afy

ENV PATH /opt/conda/bin:${PATH}
RUN pip install .