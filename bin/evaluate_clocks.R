#!/usr/bin/env Rscript

# ==============================================================================
# evaluate_clocks.R — Evaluate Ancestry Composition and Baseline Clocks
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(methylclockData)
  library(methylclock)
})

options(ExperimentHub.ask = FALSE)
options(AnnotationHub.ask = FALSE)

option_list <- list(
  make_option(c("--coefs-dir"), type = "character", default = NULL,
              help = "Directory containing model coefficient CSVs [required]", metavar = "dir", dest = "coefs_dir"),
  make_option(c("--magenta-beta"), type = "character", default = NULL,
              help = "Path to MAGENTA beta matrix RDS file [required]", metavar = "file", dest = "magenta_beta"),
  make_option(c("--magenta-meta"), type = "character", default = NULL,
              help = "Path to MAGENTA pheno metadata file (.xlsx/.csv/.rds) [required]", metavar = "file", dest = "magenta_meta"),
  make_option(c("--output-predictions"), type = "character", default = NULL,
              help = "Path to save predicted ages CSV file [required]", metavar = "file", dest = "output_predictions"),
  make_option(c("--output-summary"), type = "character", default = NULL,
              help = "Path to save summary evaluation metrics CSV file [required]", metavar = "file", dest = "output_summary")
)

parser <- OptionParser(
  usage = "%prog [options]",
  option_list = option_list,
  description = "Evaluate Baseline and Ancestry Composition DNAm age clocks on external validation cohort (e.g. MAGENTA)."
)
opt <- parse_args(parser)

# Validate required arguments
if (is.null(opt$coefs_dir) || is.null(opt$magenta_beta) || 
    is.null(opt$magenta_meta) || is.null(opt$output_predictions) || 
    is.null(opt$output_summary)) {
  print_help(parser)
  stop("Error: All arguments (--coefs-dir, --magenta-beta, --magenta-meta, --output-predictions, --output-summary) are required.", call. = FALSE)
}

cat("===========================================================\n")
cat("  Evaluating Clocks on Validation Cohort (MAGENTA)\n")
cat("===========================================================\n\n")

# Helper Functions
anti_transform_age <- function(x, adult_age = 20) {
  ifelse(x <= 0, exp(x + log(adult_age + 1)) - 1, x * (adult_age + 1) + adult_age)
}

predict_from_coefs <- function(coefs, beta_matrix, cpg_col, coef_col, intercept_val = 0, needs_anti_trafo = FALSE) {
  shared <- intersect(coefs[[cpg_col]], rownames(beta_matrix))
  if (length(shared) == 0) {
    return(rep(intercept_val, ncol(beta_matrix)))
  }
  coefs_sub <- coefs[coefs[[cpg_col]] %in% shared, , drop = FALSE]
  beta_sub <- beta_matrix[coefs_sub[[cpg_col]], , drop = FALSE]
  pred_t <- as.vector(t(coefs_sub[[coef_col]]) %*% as.matrix(beta_sub)) + intercept_val
  if (needs_anti_trafo) return(anti_transform_age(pred_t)) else return(pred_t)
}

# 1. Load Validation Data
cat("Loading validation dataset...\n")
cat("  Beta: ", opt$magenta_beta, "\n")
ext_beta <- tolower(tools::file_ext(opt$magenta_beta))
if (ext_beta == "csv") {
  beta_df <- read.csv(opt$magenta_beta, row.names = 1, check.names = FALSE)
  beta_MG <- as.matrix(beta_df)
} else {
  beta_MG <- readRDS(opt$magenta_beta)
}

cat("  Meta: ", opt$magenta_meta, "\n")
ext <- tolower(tools::file_ext(opt$magenta_meta))
if (ext %in% c("xlsx", "xls")) {
  suppressPackageStartupMessages(library(readxl))
  meta_MG <- read_excel(opt$magenta_meta)
} else if (ext == "rds") {
  meta_MG <- readRDS(opt$magenta_meta)
} else if (ext == "csv") {
  meta_MG <- read.csv(opt$magenta_meta, stringsAsFactors = FALSE)
} else {
  stop("Unsupported file extension for --magenta-meta: ", ext)
}

# Standardize / Filter metadata
if (!"Beta_ID" %in% colnames(meta_MG) && "sample_ID" %in% colnames(meta_MG)) {
  meta_MG$Beta_ID <- meta_MG$sample_ID
}
if (!"AGE_OF_EXAM" %in% colnames(meta_MG) && "age" %in% colnames(meta_MG)) {
  meta_MG$AGE_OF_EXAM <- meta_MG$age
}

meta_MG <- meta_MG %>% filter(Beta_ID %in% colnames(beta_MG), !is.na(AGE_OF_EXAM))

# Format rownames for beta
rownames(beta_MG) <- gsub("_.*", "", rownames(beta_MG))
beta_MG <- beta_MG[, meta_MG$Beta_ID]

cohort_labels <- c("CuADI" = "Cuban", "NHW" = "Non-Hisp. White", "PERUVIAN" = "Peruvian", "PRADI" = "Puerto Rican", "REAAADI" = "African American")
if ("COHORT" %in% colnames(meta_MG)) {
  meta_MG$Cohort_Label <- ifelse(meta_MG$COHORT %in% names(cohort_labels), cohort_labels[meta_MG$COHORT], meta_MG$COHORT)
} else if (!"Cohort_Label" %in% colnames(meta_MG)) {
  meta_MG$Cohort_Label <- "Unknown"
}

cat(sprintf("  Validation Cohort (All): N=%d samples x %d CpGs\n", ncol(beta_MG), nrow(beta_MG)))

datasets <- list(
  "Whole MAGENTA" = list(
    beta = beta_MG,
    age = meta_MG$AGE_OF_EXAM,
    meta = meta_MG
  )
)

if ("STATUS" %in% colnames(meta_MG) && any(meta_MG$STATUS == "CONTROL")) {
  datasets[["Controls Only"]] <- list(
    beta = beta_MG[, meta_MG$STATUS == "CONTROL"],
    age = meta_MG$AGE_OF_EXAM[meta_MG$STATUS == "CONTROL"],
    meta = meta_MG[meta_MG$STATUS == "CONTROL", ]
  )
}

# 2. Extract Baseline Clocks
cat("\nExtracting external and baseline clocks...\n")
load_DNAm_Clocks_data()

format_external <- function(coef_df, needs_anti_trafo) {
  if ("Intercept" %in% coef_df$CpGmarker) {
    intercept <- coef_df$CoefficientTraining[coef_df$CpGmarker == "Intercept"]
    cpgs <- coef_df[coef_df$CpGmarker != "Intercept", ]
  } else {
    intercept <- 0; cpgs <- coef_df
  }
  cpgs <- data.frame(cpg = cpgs$CpGmarker, coefficient = cpgs$CoefficientTraining, stringsAsFactors = FALSE)
  list(cpgs = cpgs, intercept = intercept, needs_anti_trafo = needs_anti_trafo)
}

baseline_clocks <- list(
  "Horvath" = format_external(coefHorvath, TRUE),
  "Hannum"  = format_external(coefHannum, FALSE),
  "Zhang EN"= format_external(coefEN, FALSE)
)

# Search for local baseline combined clock in --coefs-dir
local_coef_files <- c(
  file.path(opt$coefs_dir, "Combined_Full_coef.csv"),
  file.path(opt$coefs_dir, "..", "Combined_Full_coef.csv")
)
for (lcf in local_coef_files) {
  if (file.exists(lcf)) {
    local_df <- read.csv(lcf, stringsAsFactors = FALSE)
    loc_int <- local_df$coefficient[local_df$cpg == "(Intercept)"]
    loc_cpg <- local_df[local_df$cpg != "(Intercept)", ]
    baseline_clocks[["Baseline_Combined"]] <- list(cpgs = loc_cpg, intercept = loc_int, needs_anti_trafo = TRUE)
    break
  }
}

# 3. Evaluate Clocks & Store Predictions + Metrics
results_list <- list()
predictions_list <- list()

# Evaluate Baseline Clocks
for (clock_name in names(baseline_clocks)) {
  model <- baseline_clocks[[clock_name]]
  for (ds_name in names(datasets)) {
    ds <- datasets[[ds_name]]
    preds <- predict_from_coefs(model$cpgs, ds$beta, "cpg", "coefficient", model$intercept, model$needs_anti_trafo)
    
    # Store predictions
    pred_df <- data.frame(
      sample_id = colnames(ds$beta),
      actual_age = ds$age,
      cohort = ds$meta$Cohort_Label,
      dataset = ds_name,
      clock = clock_name,
      clock_type = "Baseline",
      ratio = NA_character_,
      iteration = NA_integer_,
      predicted_age = preds,
      stringsAsFactors = FALSE
    )
    predictions_list[[length(predictions_list) + 1]] <- pred_df

    # Whole dataset metric
    res_all <- data.frame(
      Dataset = ds_name, Cohort = "All", Clock_Type = "Baseline", Ratio = NA_character_, Iteration = NA_integer_, Clock = clock_name,
      MAE = mean(abs(preds - ds$age)),
      MedAE = median(abs(preds - ds$age)),
      Pearson_R = cor(preds, ds$age), stringsAsFactors = FALSE
    )
    results_list[[length(results_list) + 1]] <- res_all
    
    # Per cohort metric
    for (cohort in unique(ds$meta$Cohort_Label[!is.na(ds$meta$Cohort_Label)])) {
      idx <- ds$meta$Cohort_Label == cohort
      cohort_preds <- preds[idx]
      cohort_age <- ds$age[idx]
      if (length(cohort_preds) >= 5) {
        res_coh <- data.frame(
          Dataset = ds_name, Cohort = cohort, Clock_Type = "Baseline", Ratio = NA_character_, Iteration = NA_integer_, Clock = clock_name,
          MAE = mean(abs(cohort_preds - cohort_age)),
          MedAE = median(abs(cohort_preds - cohort_age)),
          Pearson_R = cor(cohort_preds, cohort_age), stringsAsFactors = FALSE
        )
        results_list[[length(results_list) + 1]] <- res_coh
      }
    }
  }
  cat("  Evaluated Baseline:", clock_name, "\n")
}

# Evaluate Composition Clocks
ratios <- c("100_0", "75_25", "50_50", "25_75", "0_100")
iterations <- 0:49

for (ratio in ratios) {
  cat(sprintf("\n  Evaluating Ratio %s models:\n", ratio))
  for (k in iterations) {
    # Check potential path locations
    coef_candidates <- c(
      file.path(opt$coefs_dir, sprintf("Composition_%s_iter%03d_coef.csv", ratio, k)),
      file.path(opt$coefs_dir, "composition", sprintf("Composition_%s_iter%03d_coef.csv", ratio, k))
    )
    coef_file <- coef_candidates[file.exists(coef_candidates)][1]
    if (is.na(coef_file) || !file.exists(coef_file)) {
      next
    }
    
    df <- read.csv(coef_file, stringsAsFactors = FALSE)
    intercept <- df$coefficient[df$cpg == "(Intercept)"]
    cpgs <- df[df$cpg != "(Intercept)", ]
    
    for (ds_name in names(datasets)) {
      ds <- datasets[[ds_name]]
      preds <- predict_from_coefs(cpgs, ds$beta, "cpg", "coefficient", intercept, TRUE)
      
      pred_df <- data.frame(
        sample_id = colnames(ds$beta),
        actual_age = ds$age,
        cohort = ds$meta$Cohort_Label,
        dataset = ds_name,
        clock = sprintf("Comp_%s", ratio),
        clock_type = "Composition",
        ratio = ratio,
        iteration = k,
        predicted_age = preds,
        stringsAsFactors = FALSE
      )
      predictions_list[[length(predictions_list) + 1]] <- pred_df

      # Whole dataset metric
      res_all <- data.frame(
        Dataset = ds_name, Cohort = "All", Clock_Type = "Composition", Ratio = ratio, Iteration = k, Clock = sprintf("Comp_%s", ratio),
        MAE = mean(abs(preds - ds$age)),
        MedAE = median(abs(preds - ds$age)),
        Pearson_R = cor(preds, ds$age), stringsAsFactors = FALSE
      )
      results_list[[length(results_list) + 1]] <- res_all
      
      # Per cohort metric
      for (cohort in unique(ds$meta$Cohort_Label[!is.na(ds$meta$Cohort_Label)])) {
        idx <- ds$meta$Cohort_Label == cohort
        cohort_preds <- preds[idx]
        cohort_age <- ds$age[idx]
        if (length(cohort_preds) >= 5) {
          res_coh <- data.frame(
            Dataset = ds_name, Cohort = cohort, Clock_Type = "Composition", Ratio = ratio, Iteration = k, Clock = sprintf("Comp_%s", ratio),
            MAE = mean(abs(cohort_preds - cohort_age)),
            MedAE = median(abs(cohort_preds - cohort_age)),
            Pearson_R = cor(cohort_preds, cohort_age), stringsAsFactors = FALSE
          )
          results_list[[length(results_list) + 1]] <- res_coh
        }
      }
    }
    if ((k + 1) %% 10 == 0 || k == 0) {
      cat(sprintf("    Finished iter %03d\n", k))
    }
  }
}

# 4. Save Outputs
dir.create(dirname(opt$output_predictions), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(opt$output_summary), showWarnings = FALSE, recursive = TRUE)

all_predictions <- bind_rows(predictions_list)
write.csv(all_predictions, opt$output_predictions, row.names = FALSE)
cat("\nSaved predicted ages to:", opt$output_predictions, "\n")

metrics_df <- bind_rows(results_list)
write.csv(metrics_df, opt$output_summary, row.names = FALSE)
cat("Saved summary metrics to:", opt$output_summary, "\n")

cat("\n===========================================================\n")
cat("  Evaluation Complete.\n")
cat("===========================================================\n")
