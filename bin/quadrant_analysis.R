#!/usr/bin/env Rscript

# =============================================================================
# DIVERGENCE / QUADRANT CONCORDANCE ANALYSIS
# Analyzes discordance and concordance between genetic risk (PRS) and epigenetic risk (AgeAccel).
# Performs Fisher's exact test for quadrant case/control enrichment.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
})

option_list <- list(
  make_option(c("--analysis-rds"), type = "character", default = NULL,
              help = "Path to analysis-ready dataframe RDS file", metavar = "FILE"),
  make_option(c("--output-csv"), type = "character", default = NULL,
              help = "Path to output Fisher's exact quadrant concordance stats CSV file", metavar = "FILE")
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

log_msg("=== DIVERGENCE / QUADRANT ANALYSIS ===")

# Load data
log_msg(sprintf("Loading analysis dataset from: %s", opt[["analysis-rds"]]))
df <- readRDS(opt[["analysis-rds"]])
log_msg(sprintf("Loaded dataset: %d samples", nrow(df)))

clocks <- c("EUR", "AFR", "Combined_50_50", "Horvath", "Hannum", "PhenoAge", "ZhangEN", "DunedinPACE")
clocks <- intersect(clocks, colnames(df))

prss <- intersect(c("PGS002280_std", "PGS003953_std", "PGS000823_std", "PGS003958_std"), colnames(df))
if (length(prss) == 0) {
  prss <- grep("_std$", colnames(df), value = TRUE)
}

enrichment_results <- list()

for (prs in prss) {
  df_prs <- df %>% filter(!is.na(.data[[prs]]))
  
  for (clock in clocks) {
    accel_var <- if (clock == "DunedinPACE") "DunedinPACE_centered" else paste0("AgeAccel_resid_", clock)
    if (!accel_var %in% colnames(df_prs)) next
    
    log_msg(sprintf("Running divergence analysis for PRS: %s, Clock: %s", prs, clock))
    
    # 1. Define quadrants using median split
    prs_median   <- median(df_prs[[prs]], na.rm = TRUE)
    accel_median <- median(df_prs[[accel_var]], na.rm = TRUE)
    
    df_quad <- df_prs %>%
      mutate(
        quadrant = case_when(
          .data[[prs]] >= prs_median & .data[[accel_var]] >= accel_median ~ "Q1_HighPRS_HighAccel",
          .data[[prs]] < prs_median  & .data[[accel_var]] >= accel_median ~ "Q2_LowPRS_HighAccel",
          .data[[prs]] < prs_median  & .data[[accel_var]] < accel_median  ~ "Q3_LowPRS_LowAccel",
          .data[[prs]] >= prs_median & .data[[accel_var]] < accel_median  ~ "Q4_HighPRS_LowAccel",
          TRUE ~ NA_character_
        )
      )
    
    # 2. Quadrant enrichment tests (Fisher's exact test)
    quads <- c("Q1_HighPRS_HighAccel", "Q2_LowPRS_HighAccel", "Q3_LowPRS_LowAccel", "Q4_HighPRS_LowAccel")
    for (qd in quads) {
      tbl <- table(
        factor(df_quad$quadrant == qd, levels = c("TRUE", "FALSE")),
        factor(df_quad$AD_status, levels = c("AD", "CONTROL"))
      )
      
      if (nrow(tbl) == 2 && ncol(tbl) == 2) {
        ft <- fisher.test(tbl)
        enrichment_results[[length(enrichment_results) + 1]] <- tibble(
          PRS                  = prs,
          Clock                = clock,
          quadrant             = qd,
          in_quadrant_cases    = tbl["TRUE", "AD"],
          in_quadrant_controls = tbl["TRUE", "CONTROL"],
          out_quadrant_cases   = tbl["FALSE", "AD"],
          out_quadrant_controls= tbl["FALSE", "CONTROL"],
          odds_ratio           = as.numeric(ft$estimate),
          p_value              = ft$p.value,
          conf_low             = ft$conf.int[1],
          conf_high            = ft$conf.int[2]
        )
      }
    }
  }
}

enrichment_df <- bind_rows(enrichment_results)

out_dir <- dirname(opt[["output-csv"]])
if (!dir.exists(out_dir)) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
}
write.csv(enrichment_df, opt[["output-csv"]], row.names = FALSE)
log_msg(sprintf("Completed quadrant concordance analysis. Saved %d rows to: %s",
                nrow(enrichment_df), opt[["output-csv"]]))

log_msg("=== DIVERGENCE / QUADRANT ANALYSIS COMPLETE ===")
