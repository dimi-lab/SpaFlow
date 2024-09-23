# Base image with Java, R, and Python
FROM rocker/r-ver:4.2.2

# Add the deadsnakes PPA to get Python 3.8
RUN apt-get update && apt-get install -y software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update

RUN apt-get install -y \
      pandoc \
      git \
      libcurl4-openssl-dev \
      libhdf5-serial-dev \
      libssl-dev \
      libxml2-dev \
      zlib1g-dev \
      libxt-dev \
      curl \
      apt-utils \
      openjdk-17-jdk \
      pkg-config \
      python3.8 python3.8-dev python3.8-distutils python3-pip

# Set system default python to 3.8 (10 == priority, not sure this matters much)
RUN update-alternatives --remove python /usr/bin/python2
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.8 10

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin

# Python Packages
RUN python -m pip install 'scimap==1.3.2' 'anndata==0.7.8' 'pandas==1.5.3' 'scanpy==1.9.6'

# Install R packages (ensuring dependencies are handled)
RUN R -e "install.packages(c('knitr', 'ggplot2', 'data.table', 'dplyr'))"

# FIXME: remove the repos argument
RUN R -e "install.packages('Seurat', repos = 'https://cloud.r-project.org')"
RUN R -e "install.packages(c('progressr', 'kableExtra', 'ComplexHeatmap', 'ggridges', 'clustree', 'pheatmap', 'plyr'))"
RUN R -e "install.packages(c('pander', 'CELESTA', 'Rmixmod'))" 
RUN R -e "install.packages(c('spdep', 'reshape2', 'zeallot'))"

# FIXME: move these to the apt-get above
RUN apt install vim libfontconfig1-dev -y

RUN R -e "install.packages(c('leiden', 'reticulate', 'SeuratObject'))"
RUN R -e "install.packages(c('systemfonts', 'ggforce', 'ggraph'))"
RUN R -e "install.packages(c('png'))"
# FIXME: install above
RUN R -e "install.packages('Seurat')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('ComplexHeatmap')"
RUN R -e "install.packages('clustree')"
RUN R -e "install.packages('kableExtra')"
RUN R -e "install.packages('Matrix')"


##############################################################
# CELESTA does not work–
# Seems to succeed installing, BUT, CELESTA_clustering.Rmd fails ("celesta not found")
# FIXME: move to apt above
# TODO: try with --no-install-recommends
# RUN apt-get install -y r-cran-devtools
# RUN R -e "install.packages('devtools')"
# RUN R -e "install.packages('remotes')"

# RUN R -e "remotes::install_github('plevritis-lab/CELESTA')"
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

