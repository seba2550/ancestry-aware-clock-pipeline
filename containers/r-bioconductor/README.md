# R/Bioconductor Container Environment

This directory contains the Docker configuration for the R/Bioconductor evaluation, DNA methylation clock annotation, and statistical analysis environment used in the ancestry-aware clock pipeline.

## Container Purpose
Provides a standardized R 4.4 / Bioconductor 3.19 execution environment containing specialized tools for DNA methylation analysis (`minfi`, `methylclock`), annotation packages, and statistical visualization tools.

## Included Packages

### Bioconductor Packages
- `minfi`
- `methylclock`
- `methylclockData`
- `IlluminaHumanMethylation450kanno.ilmn12.hg19`

### CRAN Packages
- `glmnet`
- `tidyverse`
- `ggplot2`
- `patchwork`
- `ggvenn`
- `ggforce`
- `pROC`
- `scales`
- `readxl`
- `optparse`
- `here`

## Build Instructions

To build the Docker image locally:

```bash
docker build -t sgonzalez/r-bioconductor:1.0.0 .
```

## Validation Instructions

To validate the container installation and ensure all required packages load correctly:

```bash
docker run --rm -v $(pwd):/app sgonzalez/r-bioconductor:1.0.0 Rscript test_container.R
```
