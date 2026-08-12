# PRS Tools Container

## Purpose
This container provides the computational environment for Polygenic Risk Score (PRS) parsing, computation, and integration with epigenetic clocks.

## Contents & Dependencies
- **Base OS**: Ubuntu 22.04 LTS
- **PLINK 1.9**: Installed at `/usr/local/bin/plink`
- **PLINK 2.0**: Installed at `/usr/local/bin/plink2`
- **R**: Base R installation
- **R Packages**:
  - `tidyverse`
  - `optparse`

## Build Command
```bash
docker build -t sgonzalez/prs-tools:1.0.0 .
```

## Validation
```bash
docker run --rm sgonzalez/prs-tools:1.0.0 /app/test_container.sh
```
