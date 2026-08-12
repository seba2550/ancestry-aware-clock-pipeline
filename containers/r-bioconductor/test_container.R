#!/usr/bin/env Rscript

# Standalone R validation script to test loaded libraries

required_packages <- c(
  "minfi",
  "methylclock",
  "methylclockData",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "glmnet",
  "tidyverse",
  "ggplot2",
  "patchwork",
  "ggvenn",
  "ggforce",
  "pROC",
  "scales",
  "readxl",
  "optparse",
  "here"
)

cat("Validating R/Bioconductor container package installation...\n\n")

failed_packages <- c()

for (pkg in required_packages) {
  cat(sprintf("Loading package: %-45s ... ", pkg))
  if (suppressPackageStartupMessages(require(pkg, character.only = TRUE))) {
    cat("OK\n")
  } else {
    cat("FAILED\n")
    failed_packages <- c(failed_packages, pkg)
  }
}

if (length(failed_packages) == 0) {
  cat("\nSUCCESS: All R/Bioconductor packages loaded successfully!\n")
  quit(status = 0)
} else {
  cat(sprintf("\nERROR: Failed to load %d package(s): %s\n", 
              length(failed_packages), paste(failed_packages, collapse = ", ")))
  quit(status = 1)
}
