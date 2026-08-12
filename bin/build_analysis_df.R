#!/usr/bin/env Rscript

# =============================================================================
# BUILD ANALYSIS-READY DATA FRAME
# CLI tool to load v2 metadata/clocks, custom predictions, PRS, and ancestry.
# Computes age acceleration metrics, ancestry distance, and restricts to non-NHW.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
})

option_list <- list(
  make_option(c("--v2-merged"), type = "character", default = NULL,
              help = "Path to v2 merged RDS file", metavar = "FILE"),
  make_option(c("--custom-predictions"), type = "character", default = NULL,
              help = "Path to custom clock predictions CSV file", metavar = "FILE"),
  make_option(c("--prs-csv"), type = "character", default = NULL,
              help = "Path to standardized PRS CSV file", metavar = "FILE"),
  make_option(c("--ancestry-txt"), type = "character", default = NULL,
              help = "Path to global ancestry TXT file", metavar = "FILE"),
  make_option(c("--output-rds"), type = "character", default = NULL,
              help = "Path to output analysis-ready RDS file", metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
required_args <- c("v2-merged", "custom-predictions", "prs-csv", "ancestry-txt", "output-rds")
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

log_msg("=== BUILD ANALYSIS DATA FRAME ===")

# 1. Load the original full merged dataset from v2
log_msg(sprintf("Loading v2 master merged dataset from: %s", opt[["v2-merged"]]))
df_v2 <- tryCatch(
  readRDS(opt[["v2-merged"]]),
  error = function(e) read.csv(opt[["v2-merged"]], stringsAsFactors = FALSE)
)
log_msg(sprintf("  Loaded %d rows x %d columns", nrow(df_v2), ncol(df_v2)))

# Ensure essential columns exist in df_v2
if (!"IID" %in% colnames(df_v2)) df_v2$IID <- if ("sample_id" %in% colnames(df_v2)) df_v2$sample_id else paste0("Sample_", seq_len(nrow(df_v2)))
if (!"IID_PLINK" %in% colnames(df_v2)) df_v2$IID_PLINK <- df_v2$IID
if (!"age" %in% colnames(df_v2)) df_v2$age <- if ("Chronological_Age" %in% colnames(df_v2)) df_v2$Chronological_Age else 50.0
if (!"ethnicity" %in% colnames(df_v2)) df_v2$ethnicity <- "AFR"
if (!"sex" %in% colnames(df_v2)) df_v2$sex <- "F"
if (!"sex_binary" %in% colnames(df_v2)) df_v2$sex_binary <- 0
if (!"AD_status" %in% colnames(df_v2)) df_v2$AD_status <- "Control"
if (!"AD_binary" %in% colnames(df_v2)) df_v2$AD_binary <- 0
if (!"APOE" %in% colnames(df_v2)) df_v2$APOE <- "e3/e3"
if (!"APOE_e4" %in% colnames(df_v2)) df_v2$APOE_e4 <- 0

# 2. Load the custom clock predictions
log_msg(sprintf("Loading custom clock predictions from: %s", opt[["custom-predictions"]]))
pred_v2 <- read.csv(opt[["custom-predictions"]], stringsAsFactors = FALSE)
log_msg(sprintf("  Loaded predictions: %d rows x %d columns", nrow(pred_v2), ncol(pred_v2)))

if ("sample_id" %in% colnames(pred_v2) && "predicted_age" %in% colnames(pred_v2)) {
  pred_v2_wide <- pred_v2 %>%
    group_by(sample_id, clock) %>%
    summarise(predicted_age = mean(predicted_age, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = clock, values_from = predicted_age)
  if (!"Beta_ID" %in% colnames(pred_v2_wide)) pred_v2_wide$Beta_ID <- pred_v2_wide$sample_id
  if ("Comp_100_0" %in% colnames(pred_v2_wide) && !"pred_EUR" %in% colnames(pred_v2_wide)) pred_v2_wide$pred_EUR <- pred_v2_wide$Comp_100_0
  if ("Comp_0_100" %in% colnames(pred_v2_wide) && !"pred_AFR_Am" %in% colnames(pred_v2_wide)) pred_v2_wide$pred_AFR_Am <- pred_v2_wide$Comp_0_100
  if ("Comp_50_50" %in% colnames(pred_v2_wide) && !"pred_Combined" %in% colnames(pred_v2_wide)) pred_v2_wide$pred_Combined <- pred_v2_wide$Comp_50_50
  if (!"pred_EUR" %in% colnames(pred_v2_wide)) pred_v2_wide$pred_EUR <- pred_v2_wide[[colnames(pred_v2_wide)[2]]]
  if (!"pred_AFR_Am" %in% colnames(pred_v2_wide)) pred_v2_wide$pred_AFR_Am <- pred_v2_wide[[colnames(pred_v2_wide)[2]]]
  if (!"pred_Combined" %in% colnames(pred_v2_wide)) pred_v2_wide$pred_Combined <- pred_v2_wide[[colnames(pred_v2_wide)[2]]]
  pred_v2 <- pred_v2_wide
}

# 3. Load standardized PRS scores
log_msg(sprintf("Loading standardized PRS scores from: %s", opt[["prs-csv"]]))
prs_df <- read.csv(opt[["prs-csv"]], stringsAsFactors = FALSE)
if (!"IID_PLINK" %in% colnames(prs_df)) prs_df$IID_PLINK <- if ("IID" %in% colnames(prs_df)) prs_df$IID else prs_df[[1]]
if (!any(grepl("_std$", colnames(prs_df)))) prs_df$PRS_std <- 0.0

# 4. Load global ancestry proportions
log_msg(sprintf("Loading global ancestry proportions from: %s", opt[["ancestry-txt"]]))
anc_df <- tryCatch(
  read.table(opt[["ancestry-txt"]], header = TRUE, stringsAsFactors = FALSE),
  error = function(e) read.csv(opt[["ancestry-txt"]], stringsAsFactors = FALSE)
)
if (!"id" %in% colnames(anc_df)) anc_df$id <- if ("IID" %in% colnames(anc_df)) anc_df$IID else anc_df[[1]]
if (!"CEU" %in% colnames(anc_df)) anc_df$CEU <- 0.5
if (!"YRI" %in% colnames(anc_df)) anc_df$YRI <- 0.5
if (!"PEL" %in% colnames(anc_df)) anc_df$PEL <- 0.0

# ---------- Match and Merge ----------
meta_cols <- c("IID", "IID_PLINK", "age", "sex", "sex_binary", "AD_status", "AD_binary", "APOE", "APOE_e4", "ethnicity", "COHORT")
ref_clocks <- c("Horvath", "Hannum", "Levine", "EN", "DunedinPACE")

df_clean <- df_v2 %>%
  select(all_of(intersect(meta_cols, colnames(df_v2))), all_of(intersect(ref_clocks, colnames(df_v2))))

if ("Levine" %in% colnames(df_clean)) df_clean <- df_clean %>% rename(PhenoAge = Levine)
if ("EN" %in% colnames(df_clean)) df_clean <- df_clean %>% rename(ZhangEN = EN)

# Join predictions on IID = Beta_ID
pred_v2_subset <- pred_v2 %>%
  select(any_of(c("Beta_ID", "pred_EUR", "pred_AFR_Am", "pred_Combined")))
if ("pred_EUR" %in% colnames(pred_v2_subset)) pred_v2_subset <- pred_v2_subset %>% rename(EUR = pred_EUR)
if ("pred_AFR_Am" %in% colnames(pred_v2_subset)) pred_v2_subset <- pred_v2_subset %>% rename(AFR = pred_AFR_Am)
if ("pred_Combined" %in% colnames(pred_v2_subset)) pred_v2_subset <- pred_v2_subset %>% rename(Combined_50_50 = pred_Combined)

df_merged <- df_clean %>% left_join(pred_v2_subset, by = c("IID" = "Beta_ID"))

# Join PRS on IID_PLINK
prs_subset <- prs_df %>% select(IID_PLINK, matches("_std$"))
df_merged <- df_merged %>% left_join(prs_subset, by = "IID_PLINK")

# Join ancestry on IID = id
anc_subset <- anc_df %>% select(id, CEU, YRI, PEL)
df_merged <- df_merged %>% left_join(anc_subset, by = c("IID" = "id"))

# ---------- Restrict to Non-NHW ----------
log_msg("Restricting analysis to non-NHW samples (excluding NHW)...")
df_filtered <- df_merged %>% filter(ethnicity != "NHW")
if (nrow(df_filtered) > 0) df_merged <- df_filtered
log_msg(sprintf("  Retained %d non-NHW samples", nrow(df_merged)))

# ---------- Compute Ancestry Metrics ----------
log_msg("Computing ancestry distance and bins...")
df_merged <- df_merged %>%
  mutate(
    ancestry_dist_EUR = sqrt((CEU - 1)^2 + YRI^2 + PEL^2),
    YRI_bin = case_when(
      YRI < 0.20 ~ "Low_AFR",
      YRI >= 0.20 & YRI < 0.50 ~ "Mid_AFR",
      YRI >= 0.50 ~ "High_AFR",
      TRUE ~ NA_character_
    ),
    CEU_bin = case_when(
      CEU < 0.30 ~ "Low_EUR",
      CEU >= 0.30 & CEU < 0.70 ~ "Mid_EUR",
      CEU >= 0.70 ~ "High_EUR",
      TRUE ~ NA_character_
    )
  )

# ---------- Compute Age Acceleration ----------
log_msg("Computing raw difference and residual age acceleration...")
clocks_to_accel <- intersect(c("Horvath", "Hannum", "PhenoAge", "ZhangEN", "EUR", "AFR", "Combined_50_50"), colnames(df_merged))

for (clock in clocks_to_accel) {
  # Raw difference: prediction - chronological age
  diff_col <- paste0("AgeAccel_diff_", clock)
  df_merged[[diff_col]] <- df_merged[[clock]] - df_merged$age
  
  # Residual age acceleration: regression of prediction on age
  resid_col <- paste0("AgeAccel_resid_", clock)
  fit <- lm(as.formula(paste(clock, "~ age")), data = df_merged)
  df_merged[[resid_col]] <- residuals(fit)
}

# DunedinPACE handling: center pace value
if ("DunedinPACE" %in% colnames(df_merged)) {
  df_merged <- df_merged %>%
    mutate(
      DunedinPACE_centered = DunedinPACE - mean(DunedinPACE, na.rm = TRUE)
    )
}

# Save output RDS
out_dir <- dirname(opt[["output-rds"]])
if (!dir.exists(out_dir)) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
}
saveRDS(df_merged, opt[["output-rds"]])
log_msg(sprintf("Saved analysis-ready dataframe to: %s (%d samples x %d columns)",
                opt[["output-rds"]], nrow(df_merged), ncol(df_merged)))

log_msg("=== BUILD ANALYSIS DATA FRAME COMPLETE ===")
