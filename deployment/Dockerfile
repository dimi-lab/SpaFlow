# Base image with Java, R, and Python
FROM rocker/r-ver:4.2.2

# Add the deadsnakes PPA to get Python 3.8
RUN apt-get update && apt-get install -y software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update

RUN apt-get install -y \
    apt-utils \
    curl \
    git \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libhdf5-serial-dev \
    libssl-dev \
    libxml2-dev \
    libxt-dev \
    openjdk-17-jdk \
    pandoc \
    pkg-config \
    python3.8 python3.8-dev python3.8-distutils python3-pip \
    vim \
    zlib1g-dev

# Set system default python to 3.8 (10 == priority, not sure this matters much)
RUN update-alternatives --remove python /usr/bin/python2
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.8 10

# Python Packages
RUN python -m pip install 'scimap==1.3.2' 'anndata==0.7.8' 'pandas==1.5.3' 'scanpy==1.9.6'

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin

# Install R packages (ensuring dependencies are handled)
RUN R --quiet -e "install.packages(c('knitr', 'ggplot2', 'data.table', 'dplyr'))"

# Notes on R dependencies:
# If we remove the 1st Seurat install, the container fails.
# If we remove the redundant installs for kableExtra, ComplexHeatmap, clustree, the container fails.
# We aren't 100% sure why this is … but … this works, so we're keeping it.
RUN R --quiet -e "install.packages('Seurat', repos = 'https://cloud.r-project.org')"
RUN R --quiet -e "install.packages(c('progressr', 'kableExtra', 'ComplexHeatmap', 'ggridges', 'clustree', 'pheatmap', 'plyr'))"
RUN R --quiet -e "install.packages(c('pander', 'CELESTA', 'Rmixmod'))"
RUN R --quiet -e "install.packages(c('spdep', 'reshape2', 'zeallot'))"
RUN R --quiet -e "install.packages(c('leiden', 'reticulate', 'SeuratObject'))"
RUN R --quiet -e "install.packages(c('systemfonts', 'ggforce', 'ggraph'))"
RUN R --quiet -e "install.packages(c('png'))"
RUN R --quiet -e "install.packages('Seurat')"
RUN R --quiet -e "install.packages('BiocManager')"
RUN R --quiet -e "BiocManager::install('ComplexHeatmap')"
RUN R --quiet -e "install.packages('clustree')"
RUN R --quiet -e "install.packages('kableExtra')"
RUN R --quiet -e "install.packages('Matrix')"


##############################################################
# CELESTA does not work–
# Seems to succeed installing, BUT, CELESTA_clustering.Rmd fails ("celesta not found")
# RUN apt-get install -y r-cran-devtools
# RUN R --quiet -e "install.packages('devtools')"
# RUN R --quiet -e "install.packages('remotes')"
# RUN R --quiet -e "remotes::install_github('plevritis-lab/CELESTA')"
##############################################################

# Cleanup to reduce image size
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

COPY data/ /app/data
COPY modules/ /app/modules
COPY scripts/ /app/scripts
COPY main.nf /app/main.nf
COPY nextflow.config /app/nextflow.config
