#!/usr/bin/env bash
set -eo pipefail

echo "=========================================="
echo "      Testing PRS Tools Container         "
echo "=========================================="

echo -n "Checking PLINK 1.9 version... "
if command -v plink &> /dev/null; then
    plink --version
else
    echo "ERROR: plink (v1.9) command not found."
    exit 1
fi

echo -n "Checking PLINK 2.0 version... "
if command -v plink2 &> /dev/null; then
    plink2 --version
else
    echo "ERROR: plink2 command not found."
    exit 1
fi

echo -n "Checking R installation... "
Rscript --version

echo "Testing R package imports (tidyverse, optparse)..."
Rscript -e "
suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})
cat('Success: R packages tidyverse and optparse loaded cleanly.\n')
"

echo "=========================================="
echo " All container tests completed successfully!"
echo "=========================================="
