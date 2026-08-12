FROM bioconductor/bioconductor_docker:3.19

# Install CRAN packages
RUN Rscript -e "install.packages(c( \
    'glmnet', \
    'tidyverse', \
    'ggplot2', \
    'patchwork', \
    'ggvenn', \
    'ggforce', \
    'pROC', \
    'scales', \
    'readxl', \
    'optparse', \
    'here' \
), repos='https://cloud.r-project.org/')"

# Install Bioconductor packages
RUN Rscript -e "BiocManager::install(c( \
    'minfi', \
    'methylclock', \
    'methylclockData', \
    'IlluminaHumanMethylation450kanno.ilmn12.hg19' \
), ask=FALSE)"

# Set working directory
WORKDIR /app

CMD ["R"]
