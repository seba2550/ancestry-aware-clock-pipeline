#!/usr/bin/env Rscript

# ==============================================================================
# characterize_cpgs.R — CpG Characterization, Annotation, and Overlap Analysis
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(minfi)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(methylclockData)
})

option_list <- list(
  make_option(c("--comp-coefs"), type = "character", default = NULL,
              help = "Path to composition coefficients CSV file [required]", metavar = "file", dest = "comp_coefs"),
  make_option(c("--learning-coefs"), type = "character", default = NULL,
              help = "Path to learning curve coefficients CSV file [required]", metavar = "file", dest = "learning_coefs"),
  make_option(c("--data-dir"), type = "character", default = NULL,
              help = "Optional path to data directory containing binary test matrices for Delta Beta calculation", metavar = "dir", dest = "data_dir"),
  make_option(c("--output-dir"), type = "character", default = NULL,
              help = "Directory to save characterization tables and RDS artifacts [required]", metavar = "dir", dest = "output_dir")
)

parser <- OptionParser(
  usage = "%prog [options]",
  option_list = option_list,
  description = "Perform comprehensive CpG characterization (genomic region, CpG island, delta beta, gene mapping, overlap stats)."
)
opt <- parse_args(parser)

if (is.null(opt$comp_coefs) || is.null(opt$learning_coefs) || is.null(opt$output_dir)) {
  print_help(parser)
  stop("Error: --comp-coefs, --learning-coefs, and --output-dir are required.", call. = FALSE)
}

cat("===========================================================\n")
cat("  Characterizing Clock CpGs\n")
cat("===========================================================\n\n")

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# Load Coefficients
cat("Loading coefficient files...\n")
comp_coefs <- read.csv(opt$comp_coefs, stringsAsFactors = FALSE)
lc_coefs   <- read.csv(opt$learning_coefs, stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 1. CpG Counts by Ratio
# ------------------------------------------------------------------------------
cat("1. Computing CpG counts by composition ratio...\n")
comp_cpgs <- comp_coefs %>%
  filter(cpg != "(Intercept)") %>%
  group_by(Ratio, Iteration) %>%
  summarise(Count = sum(coefficient != 0), .groups = "drop") %>%
  mutate(
    EUR_Proportion = case_when(
      Ratio %in% c("100/0", "100_0") ~ 1.00,
      Ratio %in% c("75/25", "75_25") ~ 0.75,
      Ratio %in% c("50/50", "50_50") ~ 0.50,
      Ratio %in% c("25/75", "25_75") ~ 0.25,
      Ratio %in% c("0/100", "0_100") ~ 0.00,
      TRUE ~ suppressWarnings(as.numeric(Ratio))
    )
  )

write.csv(comp_cpgs, file.path(opt$output_dir, "cpg_counts_by_ratio.csv"), row.names = FALSE)
cat("  Saved cpg_counts_by_ratio.csv\n")

# ------------------------------------------------------------------------------
# 2. Extract Representative Clocks & CpG Overlap Stats
# ------------------------------------------------------------------------------
cat("2. Analyzing CpG overlap across EUR, AFR, and 50/50 Clocks...\n")

if ("Clock_Type" %in% colnames(lc_coefs)) {
  eur_sub <- lc_coefs %>% filter(Clock_Type == "EUR", cpg != "(Intercept)", coefficient != 0)
  afr_sub <- lc_coefs %>% filter(Clock_Type %in% c("AFR-Am", "AFR"), cpg != "(Intercept)", coefficient != 0)
} else {
  eur_sub <- lc_coefs %>% filter(Ratio %in% c("100/0", "100_0"), cpg != "(Intercept)", coefficient != 0)
  afr_sub <- lc_coefs %>% filter(Ratio %in% c("0/100", "0_100"), cpg != "(Intercept)", coefficient != 0)
}

if ("N" %in% colnames(eur_sub)) {
  if (1650 %in% eur_sub$N) eur_sub <- eur_sub %>% filter(N == 1650) else eur_sub <- eur_sub %>% filter(N == max(N))
}
if ("Iteration" %in% colnames(eur_sub)) {
  if (29 %in% eur_sub$Iteration) eur_sub <- eur_sub %>% filter(Iteration == 29) else eur_sub <- eur_sub %>% filter(Iteration == min(Iteration))
}
cpg_eur <- unique(eur_sub$cpg)

if ("N" %in% colnames(afr_sub)) {
  if (1650 %in% afr_sub$N) afr_sub <- afr_sub %>% filter(N == 1650) else afr_sub <- afr_sub %>% filter(N == max(N))
}
if ("Iteration" %in% colnames(afr_sub)) {
  if (1 %in% afr_sub$Iteration) afr_sub <- afr_sub %>% filter(Iteration == 1) else afr_sub <- afr_sub %>% filter(Iteration == min(Iteration))
}
cpg_afr <- unique(afr_sub$cpg)

# 50/50 Combined clock
comb_sub <- comp_coefs %>% filter(Ratio %in% c("50/50", "50_50"), cpg != "(Intercept)", coefficient != 0)
if ("Iteration" %in% colnames(comb_sub)) {
  if (3 %in% comb_sub$Iteration) comb_sub <- comb_sub %>% filter(Iteration == 3) else comb_sub <- comb_sub %>% filter(Iteration == min(Iteration))
}
cpg_comb <- unique(comb_sub$cpg)

cpg_overlap_lists <- list(
  EUR = cpg_eur,
  AFR = cpg_afr,
  Combined_50_50 = cpg_comb,
  EUR_only = setdiff(cpg_eur, union(cpg_afr, cpg_comb)),
  AFR_only = setdiff(cpg_afr, union(cpg_eur, cpg_comb)),
  Combined_only = setdiff(cpg_comb, union(cpg_eur, cpg_afr)),
  EUR_AFR_shared = setdiff(intersect(cpg_eur, cpg_afr), cpg_comb),
  EUR_Combined_shared = setdiff(intersect(cpg_eur, cpg_comb), cpg_afr),
  AFR_Combined_shared = setdiff(intersect(cpg_afr, cpg_comb), cpg_eur),
  Core_3way_shared = Reduce(intersect, list(cpg_eur, cpg_afr, cpg_comb))
)

saveRDS(cpg_overlap_lists, file.path(opt$output_dir, "cpg_overlap_lists.rds"))

cpg_overlap_stats <- data.frame(
  Category = names(cpg_overlap_lists),
  Count = sapply(cpg_overlap_lists, length),
  stringsAsFactors = FALSE
)
write.csv(cpg_overlap_stats, file.path(opt$output_dir, "cpg_overlap_stats.csv"), row.names = FALSE)
cat("  Saved cpg_overlap_stats.csv and cpg_overlap_lists.rds\n")

# ------------------------------------------------------------------------------
# 3. Genomic Region & CpG Island Annotations
# ------------------------------------------------------------------------------
cat("3. Annotating CpGs with Illumina 450k genomic region & CpG island relation...\n")
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno_df <- as.data.frame(anno)[, c("Name", "UCSC_RefGene_Group", "Relation_to_Island")]
anno_df$CpG <- rownames(anno_df)

extract_primary <- function(x) {
  sapply(strsplit(as.character(x), ";"), `[`, 1)
}
anno_df$UCSC_RefGene_Group <- extract_primary(anno_df$UCSC_RefGene_Group)
anno_df$UCSC_RefGene_Group[is.na(anno_df$UCSC_RefGene_Group) | anno_df$UCSC_RefGene_Group == ""] <- "Intergenic"

anno_df <- anno_df %>%
  mutate(
    Relation_to_Island = case_when(
      Relation_to_Island == "Island" ~ "Island",
      Relation_to_Island %in% c("N_Shore", "S_Shore") ~ "Shore",
      Relation_to_Island %in% c("N_Shelf", "S_Shelf") ~ "Shelf",
      Relation_to_Island == "OpenSea" ~ "OpenSea",
      TRUE ~ as.character(Relation_to_Island)
    )
  )

horvath_df <- get_coefHorvath()
cpg_horvath <- horvath_df$CpGmarker[!grepl("Intercept", horvath_df$CpGmarker) & horvath_df$CpGmarker != ""]

hannum_df <- get_coefHannum()
cpg_hannum <- hannum_df$CpGmarker[!grepl("Intercept", hannum_df$CpGmarker) & hannum_df$CpGmarker != ""]

zhang_df <- get_coefEN()
cpg_zhang <- zhang_df$CpGmarker[!grepl("Intercept", zhang_df$CpGmarker) & zhang_df$CpGmarker != "" & grepl("^cg", zhang_df$CpGmarker)]

build_annotated_set <- function(cpgs, label) {
  df <- data.frame(CpG = cpgs, Set = label, stringsAsFactors = FALSE)
  left_join(df, anno_df, by = "CpG")
}

annotated_all <- bind_rows(
  build_annotated_set(cpg_eur, "European Clock"),
  build_annotated_set(cpg_comb, "Combined 50/50 Clock"),
  build_annotated_set(cpg_afr, "African Clock"),
  build_annotated_set(cpg_horvath, "Horvath"),
  build_annotated_set(cpg_hannum, "Hannum"),
  build_annotated_set(cpg_zhang, "Zhang EN")
)

write.csv(annotated_all, file.path(opt$output_dir, "annotated_clock_cpgs.csv"), row.names = FALSE)

# CpG Island Distribution Summary
island_dist <- annotated_all %>%
  filter(!is.na(Relation_to_Island)) %>%
  group_by(Set, Relation_to_Island) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  mutate(Proportion = Count / sum(Count))

write.csv(island_dist, file.path(opt$output_dir, "cpg_island_distribution.csv"), row.names = FALSE)

# Genomic Region Distribution Summary
region_dist <- annotated_all %>%
  filter(!is.na(UCSC_RefGene_Group)) %>%
  group_by(Set, UCSC_RefGene_Group) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  mutate(Proportion = Count / sum(Count))

write.csv(region_dist, file.path(opt$output_dir, "cpg_genomic_region_distribution.csv"), row.names = FALSE)

cat("  Saved annotated_clock_cpgs.csv, cpg_island_distribution.csv, and cpg_genomic_region_distribution.csv\n")

# ------------------------------------------------------------------------------
# 4. Gene Mapping & Core Gene List
# ------------------------------------------------------------------------------
cat("4. Performing Gene-level mapping and overlap analysis...\n")

anno_df_raw <- as.data.frame(anno)
name_col <- if ("Name" %in% colnames(anno_df_raw)) anno_df_raw$Name else rownames(anno_df_raw)
refgene_col <- if ("UCSC_RefGene_Name" %in% colnames(anno_df_raw)) as.character(anno_df_raw$UCSC_RefGene_Name) else rep("", length(name_col))

anno_gene_df <- data.frame(
  Name = name_col,
  UCSC_RefGene_Name = refgene_col,
  CpG = name_col,
  stringsAsFactors = FALSE
)

get_cpg_gene_df <- function(cpgs, clock_name) {
  df_sub <- anno_gene_df[anno_gene_df$CpG %in% cpgs & !is.na(anno_gene_df$UCSC_RefGene_Name) & anno_gene_df$UCSC_RefGene_Name != "", ]
  if (nrow(df_sub) == 0) {
    return(tibble(CpG = character(), Gene = character(), Clock = character()))
  }
  cpg_list <- character()
  gene_list <- character()
  for (i in seq_len(nrow(df_sub))) {
    genes <- unique(trimws(unlist(strsplit(as.character(df_sub$UCSC_RefGene_Name[i]), ";"))))
    genes <- genes[genes != ""]
    if (length(genes) > 0) {
      cpg_list <- c(cpg_list, rep(df_sub$CpG[i], length(genes)))
      gene_list <- c(gene_list, genes)
    }
  }
  tibble(CpG = cpg_list, Gene = gene_list, Clock = clock_name) %>% distinct()
}

eur_gene_df  <- get_cpg_gene_df(cpg_eur, "European Clock")
afr_gene_df  <- get_cpg_gene_df(cpg_afr, "African Clock")
comb_gene_df <- get_cpg_gene_df(cpg_comb, "Combined 50/50 Clock")

eur_genes  <- unique(eur_gene_df$Gene)
afr_genes  <- unique(afr_gene_df$Gene)
comb_genes <- unique(comb_gene_df$Gene)

core_genes <- Reduce(intersect, list(eur_genes, afr_genes, comb_genes))

gene_overlap_lists <- list(
  EUR_genes = eur_genes,
  AFR_genes = afr_genes,
  Combined_50_50_genes = comb_genes,
  Core_3way_shared_genes = core_genes
)
saveRDS(gene_overlap_lists, file.path(opt$output_dir, "gene_overlap_lists.rds"))

gene_overlap_stats <- data.frame(
  Category = names(gene_overlap_lists),
  Count = sapply(gene_overlap_lists, length),
  stringsAsFactors = FALSE
)
write.csv(gene_overlap_stats, file.path(opt$output_dir, "gene_overlap_stats.csv"), row.names = FALSE)

# Core gene table details
core_gene_table <- bind_rows(eur_gene_df, afr_gene_df, comb_gene_df) %>%
  filter(Gene %in% core_genes) %>%
  group_by(Gene) %>%
  summarise(
    Associated_CpGs = paste(unique(CpG), collapse = ";"),
    n_CpGs = length(unique(CpG)),
    Clocks = paste(unique(Clock), collapse = ";"),
    .groups = "drop"
  )

write.csv(core_gene_table, file.path(opt$output_dir, "core_genes.csv"), row.names = FALSE)
cat("  Saved gene_overlap_stats.csv, gene_overlap_lists.rds, and core_genes.csv\n")

# ------------------------------------------------------------------------------
# 5. Delta Beta Analysis (Optional if data-dir is specified or files exist)
# ------------------------------------------------------------------------------
if (!is.null(opt$data_dir) && dir.exists(opt$data_dir)) {
  cat("5. Calculating Delta Beta (|Delta Beta| European vs. African)...\n")
  
  load_test_means <- function(prefix, test_ids) {
    dims_file <- paste0(prefix, "_X_dims.txt")
    if (!file.exists(dims_file)) return(NULL)
    dim_lines <- readLines(dims_file)
    dim1 <- as.integer(dim_lines[1])
    dim2 <- as.integer(dim_lines[2])
    orientation <- ifelse(length(dim_lines) >= 3, dim_lines[3], "samples_cpgs")
    cpg_names <- readLines(paste0(prefix, "_cpg_names.txt"))
    cpg_names <- cpg_names[nchar(cpg_names) > 0]
    y_df <- read.csv(paste0(prefix, "_y.csv"), stringsAsFactors = FALSE)
    
    test_indices <- which(y_df$sample_id %in% test_ids)
    n_test <- length(test_indices)
    
    con <- file(paste0(prefix, "_X.bin"), "rb")
    raw_data <- readBin(con, what = "double", n = dim1 * dim2)
    close(con)
    
    if (orientation == "cpgs_samples") {
      flat_indices <- rep(test_indices, dim1) + rep((0:(dim1 - 1)) * dim2, each = n_test)
      test_data <- raw_data[flat_indices]
      rm(raw_data, flat_indices); gc(verbose = FALSE)
      test_mat <- matrix(test_data, nrow = n_test, ncol = dim1, byrow = FALSE)
    } else {
      flat_indices <- rep((test_indices - 1) * dim1, each = dim1) + rep(1:dim1, n_test)
      test_data <- raw_data[flat_indices]
      rm(raw_data, flat_indices); gc(verbose = FALSE)
      test_mat <- matrix(test_data, nrow = n_test, ncol = dim1, byrow = TRUE)
    }
    colnames(test_mat) <- cpg_names
    rownames(test_mat) <- y_df$sample_id[test_indices]
    means <- colMeans(test_mat, na.rm = TRUE)
    rm(test_mat); gc(verbose = FALSE)
    return(means)
  }

  eur_split_file <- file.path(opt$data_dir, "split_eur_test.csv")
  afr_split_file <- file.path(opt$data_dir, "split_afr_test.csv")
  
  if (file.exists(eur_split_file) && file.exists(afr_split_file)) {
    eur_test_ids <- read.csv(eur_split_file, stringsAsFactors = FALSE)$sample_id
    afr_test_ids <- read.csv(afr_split_file, stringsAsFactors = FALSE)$sample_id
    
    mean_eur <- load_test_means(file.path(opt$data_dir, "clock_combined_full_eur"), eur_test_ids)
    mean_afr <- load_test_means(file.path(opt$data_dir, "clock_combined_full_afr"), afr_test_ids)
    
    if (!is.null(mean_eur) && !is.null(mean_afr)) {
      shared_cpgs <- intersect(names(mean_eur), names(mean_afr))
      delta_beta <- abs(mean_eur[shared_cpgs] - mean_afr[shared_cpgs])
      
      df_db <- data.frame(CpG = shared_cpgs, Delta_Beta = delta_beta, stringsAsFactors = FALSE)
      write.csv(df_db, file.path(opt$output_dir, "delta_beta_per_cpg.csv"), row.names = FALSE)
      
      eur_db  <- df_db %>% filter(CpG %in% cpg_eur)  %>% mutate(Clock = "European Clock")
      afr_db  <- df_db %>% filter(CpG %in% cpg_afr)  %>% mutate(Clock = "African Clock")
      comb_db <- df_db %>% filter(CpG %in% cpg_comb) %>% mutate(Clock = "Combined 50/50 Clock")
      bg_db   <- df_db %>% mutate(Clock = "Background (All CpGs)")
      
      db_summary <- bind_rows(eur_db, afr_db, comb_db, bg_db) %>%
        group_by(Clock) %>%
        summarise(
          n_CpGs = n(),
          Mean_Delta_Beta = mean(Delta_Beta, na.rm = TRUE),
          Median_Delta_Beta = median(Delta_Beta, na.rm = TRUE),
          SD_Delta_Beta = sd(Delta_Beta, na.rm = TRUE),
          IQR_Delta_Beta = IQR(Delta_Beta, na.rm = TRUE),
          .groups = "drop"
        )
      
      write.csv(db_summary, file.path(opt$output_dir, "delta_beta_summary.csv"), row.names = FALSE)
      cat("  Saved delta_beta_per_cpg.csv and delta_beta_summary.csv\n")
    }
  }
}

cat("\n===========================================================\n")
cat("  CpG Characterization Complete.\n")
cat("  Output saved to:", opt$output_dir, "\n")
cat("===========================================================\n")
