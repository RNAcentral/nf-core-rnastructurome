FROM --platform=linux/amd64 docker.io/dincarnato/rnaframework:2.9.6

USER root

RUN apt-get update \
    && apt-get install -y --no-install-recommends r-base \
    && Rscript -e "install.packages(c('ggplot2','RColorBrewer','patchwork'), repos='https://cloud.r-project.org')" \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

USER 1000:1000
