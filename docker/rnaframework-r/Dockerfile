FROM docker.io/dincarnato/rnaframework:2.9.6

USER root

RUN apt-get update \
    && apt-get install -y --no-install-recommends r-base \
    && Rscript -e "install.packages(c('ggplot2','RColorBrewer','patchwork'), repos='https://cloud.r-project.org')" \
    && ln -sf /home/ubuntu/rnaframework/rf-count /usr/local/bin/rf-count \
    && ln -sf /home/ubuntu/rnaframework/rf-norm /usr/local/bin/rf-norm \
    && ln -sf /home/ubuntu/rnaframework/rf-fold /usr/local/bin/rf-fold \
    && ln -sf /usr/bin/R /usr/local/bin/R \
    && ln -sf /usr/bin/Rscript /usr/local/bin/Rscript \
    && ln -sf /home/ubuntu/conda/bin/RNAplot /usr/local/bin/RNAplot \
    && chmod o+rx /home/ubuntu /home/ubuntu/rnaframework /home/ubuntu/conda /home/ubuntu/conda/bin \
    && rm -rf /home/ubuntu/conda/pkgs \
    && rm -rf /home/ubuntu/conda/include \
    && rm -rf /home/ubuntu/conda/man \
    && rm -rf /home/ubuntu/conda/libexec/gcc \
    && rm -rf /home/ubuntu/conda/compiler_compat \
    && rm -rf /home/ubuntu/conda/x86_64-conda-linux-gnu \
    && rm -rf /home/ubuntu/conda/lib/python3.10 \
    && rm -rf /home/ubuntu/conda/lib/python3.12 \
    && rm -f /home/ubuntu/conda/bin/python \
    && rm -f /home/ubuntu/conda/bin/python3 \
    && rm -f /home/ubuntu/conda/bin/python3-config \
    && rm -f /home/ubuntu/conda/bin/python3.10 \
    && rm -f /home/ubuntu/conda/bin/pydoc \
    && rm -f /home/ubuntu/conda/bin/pydoc3 \
    && rm -f /home/ubuntu/conda/bin/idle3 \
    && rm -f /home/ubuntu/conda/bin/2to3 \
    && rm -f /home/ubuntu/conda/bin/bowtie \
    && rm -f /home/ubuntu/conda/bin/bowtie-align-l \
    && rm -f /home/ubuntu/conda/bin/bowtie-align-s \
    && rm -f /home/ubuntu/conda/bin/bowtie-build \
    && rm -f /home/ubuntu/conda/bin/bowtie-build-l \
    && rm -f /home/ubuntu/conda/bin/bowtie-build-s \
    && rm -f /home/ubuntu/conda/bin/bowtie-inspect \
    && rm -f /home/ubuntu/conda/bin/bowtie-inspect-l \
    && rm -f /home/ubuntu/conda/bin/bowtie-inspect-s \
    && rm -f /home/ubuntu/conda/bin/bowtie2 \
    && rm -f /home/ubuntu/conda/bin/bowtie2-align-l \
    && rm -f /home/ubuntu/conda/bin/bowtie2-align-s \
    && rm -f /home/ubuntu/conda/bin/bowtie2-build \
    && rm -f /home/ubuntu/conda/bin/bowtie2-build-l \
    && rm -f /home/ubuntu/conda/bin/bowtie2-build-s \
    && rm -f /home/ubuntu/conda/bin/bowtie2-inspect \
    && rm -f /home/ubuntu/conda/bin/bowtie2-inspect-l \
    && rm -f /home/ubuntu/conda/bin/bowtie2-inspect-s \
    && rm -f /home/ubuntu/conda/bin/bowtie2sam.pl \
    && rm -f /home/ubuntu/conda/bin/cutadapt \
    && rm -f /home/ubuntu/conda/bin/bedtools \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

USER 1000:1000
