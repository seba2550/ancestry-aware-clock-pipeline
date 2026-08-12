#!/usr/bin/env Rscript

# =============================================================================
# PRIMARY LOGISTIC REGRESSION MODELS
# Fits logistic regressions for AD ~ AgeAccel + PRS + Covariates.
# Evaluates AUC, OR, Brier Score across clocks, PRSs, and model tiers.
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
              help = "Path to output model performance metrics CSV file", metavar = "FILE")
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

log_msg("=== FIT LOGISTIC REGRESSION MODELS ===")

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
covariate_sets <- list(
  noAPOE   = c("age", "sex"),
  withAPOE = c("age", "sex", "APOE_e4")
)

results <- list()

# Helper function to compute Brier Score
brier_score <- function(observed, predicted) {
  mean((predicted - observed)^2)
}

# Helper to fit a model, extract AUC, CI, ORs
evaluate_model <- function(formula, data, clock_term = NULL, prs_term = NULL) {
  data_subset <- data %>%
    select(all_of(all.vars(formula))) %>%
    drop_na()
  
  if (nrow(data_subset) < 10) return(NULL)
  
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
  
  return(list(
    N            = nrow(data_subset),
    AUC          = auc_val,
    AUC_CI_low   = auc_ci[1],
    AUC_CI_high  = auc_ci[3],
    Brier        = brier,
    clock_OR     = clock_stats[1],
    clock_OR_low = clock_stats[2],
    clock_OR_high= clock_stats[3],
    clock_pval   = clock_stats[4],
    prs_OR       = prs_stats[1],
    prs_OR_low   = prs_stats[2],
    prs_OR_high  = prs_stats[3],
    prs_pval     = prs_stats[4]
  ))
}

log_msg("Running logistic regression model grid...")

for (prs in prss) {
  log_msg(sprintf("Processing PRS: %s", prs))
  df_prs <- df %>% filter(!is.na(.data[[prs]]))
  
  for (cov_set_name in names(covariate_sets)) {
    covs <- covariate_sets[[cov_set_name]]
    
    if ("APOE_e4" %in% covs) {
      df_prs_cov <- df_prs %>% filter(!is.na(APOE_e4))
    } else {
      df_prs_cov <- df_prs
    }
    
    for (clock in clocks) {
      methods <- if (clock == "DunedinPACE") "centered" else accel_methods
      
      for (method in methods) {
        clock_var <- if (clock == "DunedinPACE") "DunedinPACE_centered" else paste0("AgeAccel_", method, "_", clock)
        if (!clock_var %in% colnames(df_prs_cov)) next
        
        # --- Tier A: Clock-only ---
        formula_A <- as.formula(paste("AD_binary ~", paste(c(clock_var, covs), collapse = " + ")))
        res_A <- evaluate_model(formula_A, df_prs_cov, clock_term = clock_var)
        if (!is.null(res_A)) {
          results[[length(results) + 1]] <- tibble(
            clock = clock, accel_method = method, covariate_set = cov_set_name, PRS_used = prs,
            model_tier = "A_ClockOnly", N = res_A$N, AUC = res_A$AUC, AUC_CI_low = res_A$AUC_CI_low, AUC_CI_high = res_A$AUC_CI_high, Brier = res_A$Brier,
            clock_OR = res_A$clock_OR, clock_OR_low = res_A$clock_OR_low, clock_OR_high = res_A$clock_OR_high, clock_pval = res_A$clock_pval,
            prs_OR = NA, prs_OR_low = NA, prs_OR_high = NA, prs_pval = NA
          )
        }
        
        # --- Tier B: PRS-only ---
        formula_B <- as.formula(paste("AD_binary ~", paste(c(prs, covs), collapse = " + ")))
        res_B <- evaluate_model(formula_B, df_prs_cov, prs_term = prs)
        if (!is.null(res_B)) {
          results[[length(results) + 1]] <- tibble(
            clock = clock, accel_method = method, covariate_set = cov_set_name, PRS_used = prs,
            model_tier = "B_PRSOnly", N = res_B$N, AUC = res_B$AUC, AUC_CI_low = res_B$AUC_CI_low, AUC_CI_high = res_B$AUC_CI_high, Brier = res_B$Brier,
            clock_OR = NA, clock_OR_low = NA, clock_OR_high = NA, clock_pval = NA,
            prs_OR = res_B$prs_OR, prs_OR_low = res_B$prs_OR_low, prs_OR_high = res_B$prs_OR_high, prs_pval = res_B$prs_pval
          )
        }
        
        # --- Tier C: Combined ---
        formula_C <- as.formula(paste("AD_binary ~", paste(c(clock_var, prs, covs), collapse = " + ")))
        res_C <- evaluate_model(formula_C, df_prs_cov, clock_term = clock_var, prs_term = prs)
        if (!is.null(res_C)) {
          results[[length(results) + 1]] <- tibble(
            clock = clock, accel_method = method, covariate_set = cov_set_name, PRS_used = prs,
            model_tier = "C_Combined", N = res_C$N, AUC = res_C$AUC, AUC_CI_low = res_C$AUC_CI_low, AUC_CI_high = res_C$AUC_CI_high, Brier = res_C$Brier,
            clock_OR = res_C$clock_OR, clock_OR_low = res_C$clock_OR_low, clock_OR_high = res_C$clock_OR_high, clock_pval = res_C$clock_pval,
            prs_OR = res_C$prs_OR, prs_OR_low = res_C$prs_OR_low, prs_OR_high = res_C$prs_OR_high, prs_pval = res_C$prs_pval
          )
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
log_msg(sprintf("Completed logistic regressions. Saved %d results to: %s", nrow(results_df), opt[["output-csv"]]))

log_msg("=== LOGISTIC REGRESSION MODELS COMPLETE ===")
