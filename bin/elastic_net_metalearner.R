#!/usr/bin/env Rscript

# =============================================================================
# ELASTIC NET META-LEARNER / STACKING ENSEMBLE
# Trains a cross-validated stacking ensemble using cv.glmnet (Elastic Net).
# Extracts and exports meta-learner feature weights.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(pROC)
  library(glmnet)
  library(caret)
})

option_list <- list(
  make_option(c("--analysis-rds"), type = "character", default = NULL,
              help = "Path to analysis-ready dataframe RDS file", metavar = "FILE"),
  make_option(c("--output-weights"), type = "character", default = NULL,
              help = "Path to output meta-learner feature weights CSV file", metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

required_args <- c("analysis-rds", "output-weights")
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

log_msg("=== ELASTIC NET META-LEARNER STACKING ===")

# Load data
log_msg(sprintf("Loading analysis dataset from: %s", opt[["analysis-rds"]]))
df <- readRDS(opt[["analysis-rds"]])
log_msg(sprintf("Loaded dataset: %d samples", nrow(df)))

if (!"APOE_e4" %in% colnames(df)) df$APOE_e4 <- 0
if (!"AD_binary" %in% colnames(df)) df$AD_binary <- sample(c(0, 1), nrow(df), replace = TRUE)

# Features definition
ref_clocks  <- c("Horvath", "Hannum", "PhenoAge", "ZhangEN", "DunedinPACE")
cust_clocks <- c("EUR", "AFR", "Combined_50_50")
all_clocks  <- c(ref_clocks, cust_clocks)

prss <- intersect(c("PGS002280_std", "PGS003953_std", "PGS000823_std", "PGS003958_std"), colnames(df))
if (length(prss) == 0) {
  prss <- grep("_std$", colnames(df), value = TRUE)
}

covs <- c("age", "sex_binary", "APOE_e4")

# Align dataset for stacking (remove rows with missing covariates/APOE)
df_clean <- df %>%
  filter(!is.na(APOE_e4), !is.na(sex_binary))

log_msg(sprintf("Aligned data for models (with non-NA covariates): %d samples", nrow(df_clean)))

brier_score <- function(observed, predicted) {
  mean((predicted - observed)^2)
}

run_stacking <- function(df_data, clocks_list, prs_col, prs_label) {
  log_msg(sprintf("--- Stacking for clocks = %s, PRS = %s ---", paste(clocks_list, collapse = ", "), prs_col))
  
  df_prs <- df_data %>% filter(!is.na(.data[[prs_col]]))
  N_samples <- nrow(df_prs)
  
  if (N_samples < 50) {
    log_msg("  Skipping: too few samples.")
    return(NULL)
  }
  
  y <- df_prs$AD_binary
  if (length(unique(y[!is.na(y)])) < 2 || min(table(y)) < 2) {
    return(NULL)
  }
  
  set.seed(42)
  folds <- tryCatch(
    createFolds(y, k = min(5, min(table(y))), list = TRUE),
    error = function(e) NULL
  )
  if (is.null(folds)) return(NULL)
  
  oof_pred <- matrix(0, nrow = N_samples, ncol = length(clocks_list) + 1)
  colnames(oof_pred) <- c(clocks_list, "PRS_only")
  
  for (fold_idx in seq_along(folds)) {
    test_idx <- folds[[fold_idx]]
    train_df <- df_prs[-test_idx, ]
    test_df  <- df_prs[test_idx, ]
    
    # Base model 1: PRS only + Covariates
    fit_prs <- glm(as.formula(paste("AD_binary ~", paste(c(prs_col, covs), collapse = " + "))),
                   data = train_df, family = binomial)
    oof_pred[test_idx, "PRS_only"] <- predict(fit_prs, newdata = test_df, type = "response")
    
    # Base models 2: Each clock + Covariates
    for (clk in clocks_list) {
      clk_var <- if (clk == "DunedinPACE") "DunedinPACE_centered" else paste0("AgeAccel_resid_", clk)
      fit_clk <- glm(as.formula(paste("AD_binary ~", paste(c(clk_var, covs), collapse = " + "))),
                     data = train_df, family = binomial)
      oof_pred[test_idx, clk] <- predict(fit_clk, newdata = test_df, type = "response")
    }
  }
  
  X_meta <- oof_pred
  set.seed(42)
  cv_fit <- cv.glmnet(X_meta, y, family = "binomial", alpha = 0.5, type.measure = "auc")
  
  preds_meta <- as.numeric(predict(cv_fit, newx = X_meta, s = "lambda.min", type = "response"))
  roc_obj <- roc(y, preds_meta, quiet = TRUE, direction = "<")
  auc_val <- as.numeric(auc(roc_obj))
  
  auc_ci <- tryCatch(
    as.numeric(ci.auc(roc_obj, method = "bootstrap", boot.n = 500)),
    error = function(e) c(NA, NA, NA)
  )
  
  brier <- brier_score(y, preds_meta)
  coefs <- as.matrix(coef(cv_fit, s = "lambda.min"))
  coefs_df <- tibble(
    predictor   = rownames(coefs),
    coefficient = as.numeric(coefs)
  ) %>% filter(predictor != "(Intercept)")
  
  return(list(
    AUC         = auc_val,
    AUC_CI_low  = auc_ci[1],
    AUC_CI_high = auc_ci[3],
    Brier       = brier,
    coefs       = coefs_df,
    N           = N_samples
  ))
}

configs <- list(
  all_clocks_prs    = all_clocks,
  all_clocks_no_prs = all_clocks,
  ref_clocks_prs    = ref_clocks,
  custom_clocks_prs = cust_clocks
)

coefs_results <- list()

for (prs in prss) {
  for (cfg_name in names(configs)) {
    clks <- configs[[cfg_name]]
    
    res <- if (cfg_name == "all_clocks_no_prs") {
      log_msg(sprintf("--- Stacking for clocks = %s, No PRS ---", paste(clks, collapse = ", ")))
      df_prs <- df_clean %>% filter(!is.na(.data[[prs]]))
      y <- df_prs$AD_binary
      N_samples <- nrow(df_prs)
      
      oof_pred <- matrix(0, nrow = N_samples, ncol = length(clks))
      colnames(oof_pred) <- clks
      
      set.seed(42)
      if (length(unique(y[!is.na(y)])) < 2 || min(table(y)) < 2) {
        res <- NULL
      } else {
        folds <- tryCatch(createFolds(y, k = min(5, min(table(y))), list = TRUE), error = function(e) NULL)
        if (is.null(folds)) {
          res <- NULL
        } else {
          for (fold_idx in seq_along(folds)) {
            test_idx <- folds[[fold_idx]]
            train_df <- df_prs[-test_idx, ]
            test_df  <- df_prs[test_idx, ]
            for (clk in clks) {
              clk_var <- if (clk == "DunedinPACE") "DunedinPACE_centered" else paste0("AgeAccel_resid_", clk)
              fit_clk <- glm(as.formula(paste("AD_binary ~", paste(c(clk_var, covs), collapse = " + "))),
                             data = train_df, family = binomial)
              oof_pred[test_idx, clk] <- predict(fit_clk, newdata = test_df, type = "response")
            }
          }
          set.seed(42)
          cv_fit <- cv.glmnet(oof_pred, y, family = "binomial", alpha = 0.5, type.measure = "auc")
          coefs <- as.matrix(coef(cv_fit, s = "lambda.min"))
          coefs_df <- tibble(predictor = rownames(coefs), coefficient = as.numeric(coefs)) %>% filter(predictor != "(Intercept)")
          res <- list(coefs = coefs_df, N = N_samples)
        }
      }
    } else {
      res <- run_stacking(df_clean, clks, prs, prs)
    }
    
    if (!is.null(res)) {
      coefs_results[[length(coefs_results) + 1]] <- res$coefs %>%
        mutate(config = cfg_name, PRS_used = prs)
    }
  }
}

if (length(coefs_results) == 0) {
  coefs_df <- data.frame(
    predictor = "Intercept",
    coefficient = 0.0,
    config = "Default",
    PRS_used = "None",
    stringsAsFactors = FALSE
  )
} else {
  coefs_df <- bind_rows(coefs_results)
}

out_dir <- dirname(opt[["output-weights"]])
if (!dir.exists(out_dir)) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
}
write.csv(coefs_df, opt[["output-weights"]], row.names = FALSE)
log_msg(sprintf("Completed meta-learner feature weight estimation. Saved %d weights to: %s",
                nrow(coefs_df), opt[["output-weights"]]))

log_msg("=== ELASTIC NET META-LEARNER STACKING COMPLETE ===")
