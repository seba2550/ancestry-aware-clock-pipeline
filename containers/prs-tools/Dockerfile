FROM ubuntu:22.04

LABEL maintainer="sgonzalez" \
      description="PRS Tools environment container with PLINK 1.9, PLINK 2.0, R base, tidyverse, and optparse" \
      version="1.0.0"

ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies and system libraries for R packages and binary tool fetching
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    unzip \
    curl \
    ca-certificates \
    build-essential \
    git \
    gfortran \
    r-base \
    r-base-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Install PLINK 1.9
RUN wget -q https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_latest.zip -O /tmp/plink1.zip && \
    unzip /tmp/plink1.zip -d /tmp/plink1 && \
    mv /tmp/plink1/plink /usr/local/bin/plink && \
    chmod +x /usr/local/bin/plink && \
    rm -rf /tmp/plink1.zip /tmp/plink1

# Install PLINK 2.0
RUN wget -q https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip -O /tmp/plink2.zip && \
    unzip /tmp/plink2.zip -d /tmp/plink2 && \
    mv /tmp/plink2/plink2 /usr/local/bin/plink2 && \
    chmod +x /usr/local/bin/plink2 && \
    rm -rf /tmp/plink2.zip /tmp/plink2

# Install R packages: tidyverse and optparse
RUN Rscript -e "install.packages(c('optparse', 'tidyverse'), repos='https://cloud.r-project.org/', Ncpus=parallel::detectCores())"

# Set working directory
WORKDIR /app

# Copy test script into container image
COPY test_container.sh /app/test_container.sh
RUN chmod +x /app/test_container.sh

CMD ["/bin/bash"]
