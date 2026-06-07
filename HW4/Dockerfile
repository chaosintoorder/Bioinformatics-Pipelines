FROM continuumio/miniconda3:latest

RUN apt-get update && apt-get install -y \
    procps \
    && rm -rf /var/lib/apt/lists/*

COPY environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml && conda clean -a

COPY bin/plot_coverage.R /usr/local/bin/plot_coverage.R
RUN chmod +x /usr/local/bin/plot_coverage.R

ENV PATH /opt/conda/envs/hw3-pipeline/bin:$PATH

WORKDIR /data
