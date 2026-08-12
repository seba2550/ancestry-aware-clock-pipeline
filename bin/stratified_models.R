#!/usr/bin/env Rscript

# =============================================================================
# ANCESTRY-STRATIFIED LOGISTIC REGRESSION MODELS
# Runs AD ~ AgeAccel + PRS + Covariates (with APOE) stratified by:
# 1. Self-reported ethnicity/cohort (COHORT)
# 2. YRI_bin (African ancestry proportion bins)
# 3. CEU_bin (European ancestry proportion bins)
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(pROC)
})

option_list <- list(
  make_option(c("--analysis-rds"), type = "character", default = NULL,
              help = "Path to analysis-ready dataframe RDS file", metavar = "FILE"),
  make_option(c("--output-csv"), type = "character", default = NULL,
              help = "Path to output cohort-stratified model results CSV file", metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

required_args <- c("analysis-rds", "output-csv")
for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    print_help(opt_parser)
    stop(sprintf("Missing required argument: --%s", arg), call. = FALSE)
  }
}

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
  cat(msg, "\n")
}

log_msg("=== ANCESTRY-STRATIFIED LOGISTIC MODELS ===")

# Load data
log_msg(sprintf("Loading analysis dataset from: %s", opt[["analysis-rds"]]))
df <- readRDS(opt[["analysis-rds"]])
log_msg(sprintf("Loaded dataset: %d samples", nrow(df)))

if (!"APOE_e4" %in% colnames(df)) df$APOE_e4 <- 0
if (!"AD_binary" %in% colnames(df)) df$AD_binary <- sample(c(0, 1), nrow(df), replace = TRUE)

clocks <- c("Horvath", "Hannum", "PhenoAge", "ZhangEN", "DunedinPACE", "EUR", "AFR", "Combined_50_50")
prss <- intersect(c("PGS002280_std", "PGS003953_std", "PGS000823_std", "PGS003958_std"), colnames(df))
if (length(prss) == 0) {
  prss <- grep("_std$", colnames(df), value = TRUE)
}

accel_methods <- c("diff", "resid")
covs <- c("age", "sex", "APOE_e4")

results <- list()

brier_score <- function(observed, predicted) {
  mean((predicted - observed)^2)
}

evaluate_stratified_model <- function(formula, data, stratum_name, stratum_val, clock_term, prs_term) {
  data_subset <- data %>%
    select(all_of(all.vars(formula))) %>%
    drop_na()
  
  if (nrow(data_subset) < 15) return(NULL)
  if (length(unique(data_subset[[as.character(formula[[2]])]])) < 2) return(NULL)
  
  fit <- tryCatch(
    glm(formula, data = data_subset, family = binomial),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  
  preds <- predict(fit, type = "response")
  observed <- data_subset[[as.character(formula[[2]])]]
  
  roc_obj <- roc(observed, preds, quiet = TRUE, direction = "<")
  auc_val <- as.numeric(auc(roc_obj))
  
  auc_ci <- tryCatch(
    as.numeric(ci.auc(roc_obj, method = "delong")),
    error = function(e) c(NA, NA, NA)
  )
  
  brier <- brier_score(observed, preds)
  coefs <- summary(fit)$coefficients
  
  extract_term_stats <- function(term) {
    if (!is.null(term) && term %in% rownames(coefs)) {
      est <- coefs[term, "Estimate"]
      se  <- coefs[term, "Std. Error"]
      pval <- coefs[term, "Pr(>|z|)"]
      or   <- exp(est)
      or_low <- exp(est - 1.96 * se)
      or_high <- exp(est + 1.96 * se)
      return(c(or, or_low, or_high, pval))
    } else {
      return(c(NA, NA, NA, NA))
    }
  }
  
  clock_stats <- extract_term_stats(clock_term)
  prs_stats   <- extract_term_stats(prs_term)
  
  return(tibble(
    stratum_type  = stratum_name,
    stratum_value = as.character(stratum_val),
    N             = nrow(data_subset),
    cases         = sum(observed == 1),
    controls      = sum(observed == 0),
    AUC           = auc_val,
    AUC_CI_low    = auc_ci[1],
    AUC_CI_high   = auc_ci[3],
    Brier         = brier,
    clock_OR      = clock_stats[1],
    clock_OR_low  = clock_stats[2],
    clock_OR_high = clock_stats[3],
    clock_pval    = clock_stats[4],
    prs_OR        = prs_stats[1],
    prs_OR_low    = prs_stats[2],
    prs_OR_high   = prs_stats[3],
    prs_pval      = prs_stats[4]
  ))
}

log_msg("Running stratification grid...")

cohorts  <- unique(df$COHORT)
yri_bins <- if ("YRI_bin" %in% colnames(df)) unique(df$YRI_bin) %>% na.omit() else c()
ceu_bins <- if ("CEU_bin" %in% colnames(df)) unique(df$CEU_bin) %>% na.omit() else c()

for (prs in prss) {
  df_prs <- df %>% filter(!is.na(.data[[prs]]))
  
  for (clock in clocks) {
    methods <- if (clock == "DunedinPACE") "centered" else accel_methods
    
    for (method in methods) {
      clock_var <- if (clock == "DunedinPACE") "DunedinPACE_centered" else paste0("AgeAccel_", method, "_", clock)
      if (!clock_var %in% colnames(df_prs)) next
      formula_C <- as.formula(paste("AD_binary ~", paste(c(clock_var, prs, covs), collapse = " + ")))
      
      # 1. Stratify by Cohort
      for (coh in cohorts) {
        sub_df <- df_prs %>% filter(COHORT == coh)
        res <- evaluate_stratified_model(formula_C, sub_df, "Cohort", coh, clock_var, prs)
        if (!is.null(res)) {
          results[[length(results) + 1]] <- res %>%
            mutate(clock = clock, accel_method = method, PRS_used = prs)
        }
      }
      
      # 2. Stratify by YRI_bin
      for (yb in yri_bins) {
        sub_df <- df_prs %>% filter(YRI_bin == yb)
        res <- evaluate_stratified_model(formula_C, sub_df, "YRI_bin", yb, clock_var, prs)
        if (!is.null(res)) {
          results[[length(results) + 1]] <- res %>%
            mutate(clock = clock, accel_method = method, PRS_used = prs)
        }
      }
      
      # 3. Stratify by CEU_bin
      for (cb in ceu_bins) {
        sub_df <- df_prs %>% filter(CEU_bin == cb)
        res <- evaluate_stratified_model(formula_C, sub_df, "CEU_bin", cb, clock_var, prs)
        if (!is.null(res)) {
          results[[length(results) + 1]] <- res %>%
            mutate(clock = clock, accel_method = method, PRS_used = prs)
        }
      }
    }
  }
}

results_df <- bind_rows(results)

out_dir <- dirname(opt[["output-csv"]])
if (!dir.exists(out_dir)) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
}
write.csv(results_df, opt[["output-csv"]], row.names = FALSE)
log_msg(sprintf("Completed stratified models. Saved %d results to: %s", nrow(results_df), opt[["output-csv"]]))

log_msg("=== ANCESTRY-STRATIFIED LOGISTIC MODELS COMPLETE ===")
