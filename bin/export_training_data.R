#!/usr/bin/env Rscript

# ==============================================================================
# export_training_data.R — Export full pooled EUR+AFR data for Python training
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("--eur-beta"), type = "character", default = NULL,
              help = "Path to EUR beta matrix RDS file [required]", metavar = "file", dest = "eur_beta"),
  make_option(c("--eur-meta"), type = "character", default = NULL,
              help = "Path to EUR metadata RDS file [required]", metavar = "file", dest = "eur_meta"),
  make_option(c("--afr-beta"), type = "character", default = NULL,
              help = "Path to AFR beta matrix RDS file [required]", metavar = "file", dest = "afr_beta"),
  make_option(c("--afr-meta"), type = "character", default = NULL,
              help = "Path to AFR metadata RDS file [required]", metavar = "file", dest = "afr_meta"),
  make_option(c("--output-prefix"), type = "character", default = NULL,
              help = "Path prefix for binary matrix export (e.g. /path/to/clock_combined_full) [required]", metavar = "prefix", dest = "output_prefix")
)

parser <- OptionParser(
  usage = "%prog [options]",
  option_list = option_list,
  description = "Export pooled European and African DNA methylation datasets to binary matrix format for downstream Python training."
)
opt <- parse_args(parser)

# Validate required arguments
if (is.null(opt$eur_beta) || is.null(opt$eur_meta) || 
    is.null(opt$afr_beta) || is.null(opt$afr_meta) || 
    is.null(opt$output_prefix)) {
  print_help(parser)
  stop("Error: --eur-beta, --eur-meta, --afr-beta, --afr-meta, and --output-prefix are required.", call. = FALSE)
}

cat("===========================================================\n")
cat("  Exporting FULL Combined (EUR+AFR) Data (No Splits)\n")
cat("===========================================================\n\n")

# Age transformation formula
transform_age <- function(age, adult_age = 20) {
  ifelse(age <= adult_age, log(age + 1) - log(adult_age + 1), (age - adult_age) / (adult_age + 1))
}

# Ensure destination directory exists
prefix <- opt$output_prefix
out_dir <- dirname(prefix)
if (out_dir != "." && out_dir != "") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
}

load_object <- function(file_path) {
  if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
    df <- read.csv(file_path, check.names = FALSE)
    if ("CpG" %in% colnames(df)) {
      mat <- as.matrix(df[, -1, drop = FALSE])
      rownames(mat) <- df[[1]]
      return(mat)
    } else {
      mat <- as.matrix(df)
      return(mat)
    }
  } else {
    return(readRDS(file_path))
  }
}

load_meta <- function(file_path) {
  if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
    df <- read.csv(file_path, stringsAsFactors = FALSE)
    if (!"sample_ID" %in% colnames(df) && "sample_id" %in% colnames(df)) {
      df$sample_ID <- df$sample_id
    }
    return(df)
  } else {
    return(readRDS(file_path))
  }
}

# 1. Load EUR
cat("Loading EUR dataset from:\n  Beta: ", opt$eur_beta, "\n  Meta: ", opt$eur_meta, "\n")
beta_A <- load_object(opt$eur_beta)
meta_A <- load_meta(opt$eur_meta)
meta_A <- meta_A[match(colnames(beta_A), meta_A$sample_ID), ]

# Apply explicit EUR age filter >= 30 if sufficient samples exist
if (sum(meta_A$age >= 30) > 0) {
  keep_eur <- meta_A$age >= 30
  meta_A <- meta_A[keep_eur, ]
  beta_A <- beta_A[, keep_eur, drop = FALSE]
}

cat(sprintf("  EUR (Filtered >= 30y): %d samples x %d CpGs\n", ncol(beta_A), nrow(beta_A)))

# 2. Load AFR
cat("Loading AFR dataset from:\n  Beta: ", opt$afr_beta, "\n  Meta: ", opt$afr_meta, "\n")
beta_B <- load_object(opt$afr_beta)
meta_B <- load_meta(opt$afr_meta)
meta_B <- meta_B[match(colnames(beta_B), meta_B$sample_ID), ]
cat(sprintf("  AFR: %d samples x %d CpGs\n", ncol(beta_B), nrow(beta_B)))

# 3. Combine
shared_cpgs <- intersect(rownames(beta_A), rownames(beta_B))
cat(sprintf("  Shared CpGs: %d\n", length(shared_cpgs)))

# Memory mindful load
beta_combined <- cbind(beta_A[shared_cpgs, ], beta_B[shared_cpgs, ])
meta_combined <- rbind(
  data.frame(sample_ID = meta_A$sample_ID, age = meta_A$age, stringsAsFactors = FALSE),
  data.frame(sample_ID = meta_B$sample_ID, age = meta_B$age, stringsAsFactors = FALSE)
)

# Free individual matrices immediately to save max overhead
rm(beta_A, beta_B, meta_A, meta_B)
gc(verbose = FALSE)

n_samples <- ncol(beta_combined)
n_features <- nrow(beta_combined)
cat(sprintf("\n  Combined Pool: %d samples x %d CpGs\n", n_samples, n_features))

# 4. Transpose
cat("  Transposing matrix for Python...\n")
X_train <- t(beta_combined)
rm(beta_combined)
gc(verbose = FALSE)

y_raw <- meta_combined$age
y_transformed <- transform_age(y_raw)

# 5. Export
cat("  Writing binary X.bin (Sequential double vectors)...\n")
con <- file(paste0(prefix, "_X.bin"), "wb")
for (i in seq_len(n_samples)) {
  writeBin(as.double(X_train[i, ]), con)
}
close(con)
cat("  Binary export complete!\n")

writeLines(c(as.character(n_samples), as.character(n_features)), paste0(prefix, "_X_dims.txt"))
writeLines(colnames(X_train), paste0(prefix, "_cpg_names.txt"))
write.csv(
  data.frame(sample_id = meta_combined$sample_ID, age_raw = y_raw, age_transformed = y_transformed),
  paste0(prefix, "_y.csv"),
  row.names = FALSE
)

file_size_gb <- file.info(paste0(prefix, "_X.bin"))$size / 1e9
cat(sprintf("  Binary file size: %.2f GB\n", file_size_gb))
cat(sprintf("  Exported to: %s_*\n\n", prefix))

cat("===========================================================\n")
cat("  Memory-safe export complete.\n")
cat("===========================================================\n")
