# Base image with Java, R, and Python
FROM rocker/tidyverse:4.2.2

# Install system dependencies
RUN apt-get update && apt-get install -y \
    pandoc \
    git \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    libxt-dev \
    curl \
    apt-utils \
    software-properties-common

# Add the deadsnakes PPA to get Python 3.8
RUN add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update

# Java Installation (OpenJDK 17)
RUN apt-get install -y openjdk-17-jdk

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin

# Install Python 
RUN apt-get install -y python3.8 python3-pip

# Python Packages
RUN pip3 install 'scimap==1.3.2' 'anndata==0.7.8' 'pandas==1.5.3' 'scanpy==1.9.6'    

# Install R packages (ensuring dependencies are handled)
RUN R -e "install.packages(c('knitr', 'ggplot2', 'data.table', 'dplyr'))"
RUN R -e "install.packages('Seurat', repos = 'https://cloud.r-project.org')"
RUN R -e "install.packages(c('progressr', 'kableExtra', 'ComplexHeatmap', 'ggridges', 'clustree', 'pheatmap', 'plyr'))"
RUN R -e "install.packages(c('pander', 'CELESTA', 'Rmixmod'))" 
RUN R -e "install.packages(c('spdep', 'reshape2', 'zeallot'))"

# Cleanup to reduce image size
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

#Expose the default port
EXPOSE 8080
