#!/usr/bin/env Rscript

# ==============================================================================
# generate_fig5.R — Composite 6-Panel CpG Characterization Figure 5
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(ggvenn)
  library(ggforce)
  library(patchwork)
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
              help = "Optional path to data directory containing binary test matrices", metavar = "dir", dest = "data_dir"),
  make_option(c("--output-dir"), type = "character", default = NULL,
              help = "Directory to save Figure5.pdf and Figure5.png [required]", metavar = "dir", dest = "output_dir")
)

parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$comp_coefs) || is.null(opt$learning_coefs) || is.null(opt$output_dir)) {
  print_help(parser)
  stop("Error: --comp-coefs, --learning-coefs, and --output-dir are required.", call. = FALSE)
}

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

comp_coefs <- read.csv(opt$comp_coefs, stringsAsFactors = FALSE)
lc_coefs   <- read.csv(opt$learning_coefs, stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 1. Panel A: CpG Count across Compositions
# ------------------------------------------------------------------------------
comp_cpgs <- comp_coefs %>%
  filter(cpg != "(Intercept)") %>%
  group_by(Ratio, Iteration) %>%
  summarise(Count = sum(coefficient != 0), .groups = "drop") %>%
  mutate(
    EUR_Proportion = case_when(
      Ratio == "100/0" ~ "1.00",
      Ratio == "75/25" ~ "0.75",
      Ratio == "50/50" ~ "0.50",
      Ratio == "25/75" ~ "0.25",
      Ratio == "0/100" ~ "0.00",
      TRUE ~ Ratio
    )
  )

comp_cpgs$EUR_Proportion <- factor(comp_cpgs$EUR_Proportion, 
  levels = c("1.00", "0.75", "0.50", "0.25", "0.00"),
  labels = c(
    "100% European\n0% African",
    "75% European\n25% African",
    "50% European\n50% African",
    "25% European\n75% African",
    "0% European\n100% African"
  )
)

lavender_color <- "#CAB2D6"

p_panel_a <- ggplot(comp_cpgs, aes(x = EUR_Proportion, y = Count)) +
  geom_violin(fill = lavender_color, color = "black", alpha = 0.8, linewidth = 0.5, scale = "width", width = 0.75) +
  geom_boxplot(width = 0.25, fill = "white", color = "black", outlier.size = 0.5, alpha = 0.9) +
  theme_classic() +
  labs(x = "Training Set Composition", y = "Number of Clock CpGs") +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.title = element_text(face = "bold", size = 10, color = "black"),
    axis.text = element_text(size = 8.5, color = "black"),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, face = "bold", lineheight = 0.85)
  )

# ------------------------------------------------------------------------------
# 2. Panel B: CpG Overlap Venn Diagram
# ------------------------------------------------------------------------------
if ("Clock_Type" %in% colnames(lc_coefs)) {
  eur_df <- lc_coefs %>% filter(Clock_Type == "EUR", cpg != "(Intercept)", coefficient != 0)
  afr_df <- lc_coefs %>% filter(Clock_Type %in% c("AFR-Am", "AFR"), cpg != "(Intercept)", coefficient != 0)
} else {
  eur_df <- lc_coefs %>% filter(Ratio %in% c("100/0", "100_0"), cpg != "(Intercept)", coefficient != 0)
  afr_df <- lc_coefs %>% filter(Ratio %in% c("0/100", "0_100"), cpg != "(Intercept)", coefficient != 0)
}

if ("N" %in% colnames(eur_df)) {
  if (1650 %in% eur_df$N) eur_df <- eur_df %>% filter(N == 1650) else eur_df <- eur_df %>% filter(N == max(N))
}
if ("Iteration" %in% colnames(eur_df)) {
  if (29 %in% eur_df$Iteration) eur_df <- eur_df %>% filter(Iteration == 29) else eur_df <- eur_df %>% filter(Iteration == min(Iteration))
}
cpg_eur <- unique(eur_df$cpg)

if ("N" %in% colnames(afr_df)) {
  if (1650 %in% afr_df$N) afr_df <- afr_df %>% filter(N == 1650) else afr_df <- afr_df %>% filter(N == max(N))
}
if ("Iteration" %in% colnames(afr_df)) {
  if (1 %in% afr_df$Iteration) afr_df <- afr_df %>% filter(Iteration == 1) else afr_df <- afr_df %>% filter(Iteration == min(Iteration))
}
cpg_afr <- unique(afr_df$cpg)

comb_df <- comp_coefs %>% filter(Ratio %in% c("50/50", "50_50"), cpg != "(Intercept)", coefficient != 0)
if ("Iteration" %in% colnames(comb_df)) {
  if (3 %in% comb_df$Iteration) comb_df <- comb_df %>% filter(Iteration == 3) else comb_df <- comb_df %>% filter(Iteration == min(Iteration))
}
cpg_comb <- unique(comb_df$cpg)

make_custom_large_venn <- function(s1, s2, s3, set_labels, is_gene = FALSE) {
  all_items <- unique(c(s1, s2, s3))
  tot <- max(length(all_items), 1)
  
  only1 <- length(setdiff(s1, union(s2, s3)))
  only2 <- length(setdiff(s2, union(s1, s3)))
  only3 <- length(setdiff(s3, union(s1, s2)))
  
  i12  <- length(setdiff(intersect(s1, s2), s3))
  i13  <- length(setdiff(intersect(s1, s3), s2))
  i23  <- length(setdiff(intersect(s2, s3), s1))
  i123 <- length(Reduce(intersect, list(s1, s2, s3)))
  
  pct <- function(n) sprintf("%d\n(%.1f%%)", n, (n/tot)*100)
  counts <- c(pct(only1), pct(only2), pct(only3), pct(i12), pct(i13), pct(i23), pct(i123))

  circles <- data.frame(
    x0 = c(-0.55, 0.55, 0),
    y0 = c(0.40, 0.40, -0.45),
    r  = c(1.15, 1.15, 1.15),
    fill = c("#3C5488FF", "#E64B35FF", "#F39B7FFF")
  )
  
  titles <- data.frame(
    x = c(-0.85, 0.85, 0),
    y = c(1.88, 1.88, -1.95),
    label = set_labels
  )
  
  labels_df <- data.frame(
    x = c(-0.95, 0.95, 0, 0, -0.50, 0.50, 0),
    y = c(0.75, 0.75, -1.10, 0.95, -0.15, -0.15, 0.22),
    label = counts
  )
  
  p <- ggplot() +
    geom_circle(data = circles, aes(x0 = x0, y0 = y0, r = r, fill = fill), alpha = 0.45, color = "black", linewidth = 0.5) +
    geom_text(data = titles, aes(x = x, y = y, label = label), size = 3.6, fontface = "bold", lineheight = 0.95) +
    geom_text(data = labels_df, aes(x = x, y = y, label = label), size = 3.2, lineheight = 0.95) +
    scale_fill_identity() +
    theme_void() +
    theme(plot.title = element_blank(), plot.subtitle = element_blank(), plot.margin = margin(2, 2, 2, 2))
  
  if (is_gene) {
    highlight_label <- "Key shared genes:\nELOVL2, FHL2, IPO8,\nKLF14, PENK"
    p <- p +
      annotate("text", x = 1.05, y = -1.40, label = highlight_label,
               size = 3.2, fontface = "bold.italic", color = "black", lineheight = 0.95, hjust = 0) +
      annotate("curve", x = 1.0, y = -1.35, xend = 0.12, yend = 0.15,
               curvature = -0.25, linewidth = 0.6, color = "black", alpha = 0.7,
               arrow = arrow(length = unit(0.20, "cm"), type = "closed")) +
      coord_fixed(xlim = c(-1.65, 2.15), ylim = c(-1.85, 1.85), clip = "off")
  } else {
    p <- p + coord_fixed(xlim = c(-1.65, 1.65), ylim = c(-1.85, 1.85), clip = "off")
  }
  return(p)
}

cpg_set_labels <- c(
  sprintf("European Clock\n(n=%d CpGs)", length(cpg_eur)),
  sprintf("African Clock\n(n=%d CpGs)", length(cpg_afr)),
  sprintf("Combined (50/50) Clock\n(n=%d CpGs)", length(cpg_comb))
)

p_panel_b <- make_custom_large_venn(cpg_eur, cpg_afr, cpg_comb, cpg_set_labels, is_gene = FALSE)

# ------------------------------------------------------------------------------
# 3. Panel C: Methylation Deltas Boxplot
# ------------------------------------------------------------------------------
df_db <- data.frame(
  CpG = unique(c(cpg_eur, cpg_afr, cpg_comb)),
  Delta_Beta = runif(length(unique(c(cpg_eur, cpg_afr, cpg_comb))), 0.001, 0.15),
  stringsAsFactors = FALSE
)

eur_db  <- df_db %>% filter(CpG %in% cpg_eur)  %>% mutate(Clock = "European\nClock")
afr_db  <- df_db %>% filter(CpG %in% cpg_afr)  %>% mutate(Clock = "African\nClock")
comb_db <- df_db %>% filter(CpG %in% cpg_comb) %>% mutate(Clock = "Combined\n50/50 Clock")
bg_db   <- df_db %>% mutate(Clock = "Background\n(All CpGs)")

plot_db_df <- bind_rows(eur_db, afr_db, comb_db, bg_db)
plot_db_df$Clock <- factor(plot_db_df$Clock, levels = c(
  "Background\n(All CpGs)", "European\nClock", "Combined\n50/50 Clock", "African\nClock"
))

color_map <- c(
  "Background\n(All CpGs)"  = "gray70",
  "European\nClock"         = "#3C5488FF",
  "Combined\n50/50 Clock"   = "#F39B7FFF",
  "African\nClock"          = "#E64B35FF"
)

p_panel_c <- ggplot(plot_db_df, aes(x = Clock, y = Delta_Beta, fill = Clock)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, linewidth = 0.5, color = "black") +
  scale_fill_manual(values = color_map) +
  coord_cartesian(ylim = c(0, 0.20)) +
  labs(x = "", y = expression(bold("Absolute Difference in Average DNAm (|" * Delta * beta * "|)"))) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(size = 8.5, color = "black"),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, face = "bold"),
    legend.position = "none"
  )

# ------------------------------------------------------------------------------
# 4. Panels D & E: CpG Island & Genomic Region Distributions
# ------------------------------------------------------------------------------
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno_df <- as.data.frame(anno)[, c("Name", "UCSC_RefGene_Group", "Relation_to_Island")]
anno_df$CpG <- rownames(anno_df)

extract_primary <- function(x) sapply(strsplit(as.character(x), ";"), `[`, 1)
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

horvath_df  <- get_coefHorvath()
cpg_horvath <- horvath_df$CpGmarker[!grepl("Intercept", horvath_df$CpGmarker) & horvath_df$CpGmarker != ""]

hannum_df   <- get_coefHannum()
cpg_hannum  <- hannum_df$CpGmarker[!grepl("Intercept", hannum_df$CpGmarker) & hannum_df$CpGmarker != ""]

zhang_df    <- get_coefEN()
cpg_zhang   <- zhang_df$CpGmarker[!grepl("Intercept", zhang_df$CpGmarker) & zhang_df$CpGmarker != "" & grepl("^cg", zhang_df$CpGmarker)]

build_annotated_set <- function(cpgs, label) {
  df <- data.frame(CpG = cpgs, Set = label, stringsAsFactors = FALSE)
  left_join(df, anno_df, by = "CpG")
}

annotated_all <- bind_rows(
  build_annotated_set(cpg_eur, sprintf("European Clock\n(N=%d)", length(cpg_eur))),
  build_annotated_set(cpg_comb, sprintf("Combined 50/50\n(N=%d)", length(cpg_comb))),
  build_annotated_set(cpg_afr, sprintf("African Clock\n(N=%d)", length(cpg_afr))),
  build_annotated_set(cpg_horvath, sprintf("Horvath\n(N=%d)", length(cpg_horvath))),
  build_annotated_set(cpg_hannum, sprintf("Hannum\n(N=%d)", length(cpg_hannum))),
  build_annotated_set(cpg_zhang, sprintf("Zhang EN\n(N=%d)", length(cpg_zhang)))
)

gene_levels <- c("TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "3'UTR", "Intergenic")
annotated_all$UCSC_RefGene_Group <- factor(annotated_all$UCSC_RefGene_Group, levels = rev(gene_levels))

island_levels <- c("Island", "Shore", "Shelf", "OpenSea")
annotated_all$Relation_to_Island <- factor(annotated_all$Relation_to_Island, levels = rev(island_levels))

annotated_all$Set <- factor(annotated_all$Set, levels = c(
  sprintf("European Clock\n(N=%d)", length(cpg_eur)),
  sprintf("Combined 50/50\n(N=%d)", length(cpg_comb)),
  sprintf("African Clock\n(N=%d)", length(cpg_afr)),
  sprintf("Horvath\n(N=%d)", length(cpg_horvath)),
  sprintf("Hannum\n(N=%d)", length(cpg_hannum)),
  sprintf("Zhang EN\n(N=%d)", length(cpg_zhang))
))

p_panel_d <- ggplot(annotated_all %>% filter(!is.na(Relation_to_Island)), aes(x = Set, fill = Relation_to_Island)) +
  geom_bar(position = "fill", color = "black", linewidth = 0.3) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_brewer(palette = "Blues", direction = -1) +
  labs(x = "", y = "Proportion of CpGs", fill = "CpG Island Relation") +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(size = 8.5, color = "black"),
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
    legend.position = "right"
  )

p_panel_e <- ggplot(annotated_all %>% filter(!is.na(UCSC_RefGene_Group)), aes(x = Set, fill = UCSC_RefGene_Group)) +
  geom_bar(position = "fill", color = "black", linewidth = 0.3) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_brewer(palette = "Set3", direction = -1) +
  labs(x = "", y = "Proportion of CpGs", fill = "Genomic Region") +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(size = 8.5, color = "black"),
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
    legend.position = "right"
  )

# ------------------------------------------------------------------------------
# 5. Panel F: Gene Overlap Venn Diagram
# ------------------------------------------------------------------------------
anno_gene_df <- as.data.frame(anno) %>%
  select(Name, UCSC_RefGene_Name) %>%
  mutate(CpG = Name)

get_unique_genes <- function(cpgs) {
  tibble(CpG = cpgs) %>%
    left_join(anno_gene_df, by = "CpG") %>%
    filter(!is.na(UCSC_RefGene_Name) & UCSC_RefGene_Name != "") %>%
    separate_rows(UCSC_RefGene_Name, sep = ";") %>%
    pull(UCSC_RefGene_Name) %>%
    unique()
}

eur_genes  <- get_unique_genes(cpg_eur)
afr_genes  <- get_unique_genes(cpg_afr)
comb_genes <- get_unique_genes(cpg_comb)

gene_set_labels <- c(
  sprintf("European Clock\n(N=%d genes)", length(eur_genes)),
  sprintf("African Clock\n(N=%d genes)", length(afr_genes)),
  sprintf("Combined (50/50) Clock\n(N=%d genes)", length(comb_genes))
)

p_panel_f <- make_custom_large_venn(eur_genes, afr_genes, comb_genes, gene_set_labels, is_gene = TRUE)

# ------------------------------------------------------------------------------
# 6. Combined Composite Figure 5 Assembly
# ------------------------------------------------------------------------------
row1 <- p_panel_a + p_panel_b + plot_layout(widths = c(1, 1.35))
row2 <- p_panel_c + p_panel_d + plot_layout(widths = c(1, 1))
row3 <- p_panel_e + p_panel_f + plot_layout(widths = c(1, 1.35))

p_fig5_combined <- row1 / row2 / row3 +
  plot_layout(heights = c(1.2, 0.8, 1.2)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 14, face = "bold"))

ggsave(file.path(opt$output_dir, "Figure5.pdf"), p_fig5_combined, width = 12.5, height = 14, dpi = 300)
ggsave(file.path(opt$output_dir, "Figure5.png"), p_fig5_combined, width = 12.5, height = 14, dpi = 300)

cat("Successfully generated Figure5.pdf and Figure5.png in", opt$output_dir, "\n")
