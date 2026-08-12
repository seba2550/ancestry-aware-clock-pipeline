# Pipeline Output Directory Structure & Output Descriptions

All published files produced by the **Ancestry-Aware Methylation Clock Nextflow Pipeline** are saved under the output directory specified by `--outdir` (default: `./results`).

---

## Directory Overview

```
results/
├── data_prep/
│   └── clock_combined_full*
├── models/
│   ├── composition_v2_coefs.csv
│   ├── Composition_summary.csv
│   ├── full_clock_coef.csv
│   └── full_clock_summary.csv
├── characterization/
│   ├── cpg_characterization_*
│   └── core_genes_custom.csv
├── evaluation/
│   ├── magenta_all_predicted_ages_and_acceleration.csv
│   └── composition_evaluation_summary.csv
├── figures/
│   ├── Figure5.pdf
│   └── Figure5.png
└── aim3_models/
    ├── primary_logistic_models.csv
    ├── stratified_logistic_models.csv
    ├── metalearner_feature_weights.csv
    └── quadrant_concordance_stats.csv
```

---

## Output Modules & Detailed Descriptions

### 1. Data Preparation (`results/data_prep/`)

**Process**: `EXPORT_TRAINING_DATA`

Contains unified binary methylation matrices and sample metadata files created by merging European (EUR, GSE55763) and African (AFR, GSE210255) training cohorts.

| File Pattern | Format | Description |
|---|---|---|
| `clock_combined_full*` | Binary Matrix / RDS | Harmonized training methylation beta-value matrix and phenotype metadata used as input for model training. |

---

### 2. Clock Models (`results/models/`)

**Processes**: `COLLECT_TRAINING_RESULTS`, `TRAIN_FULL_CLOCK`

Contains trained Elastic Net model coefficients and summary metrics for both composition gradient models (100:0, 75:25, 50:50, 25:75, 0:100 EUR:AFR) and the full unified clock model.

| File | Format | Description |
|---|---|---|
| `composition_v2_coefs.csv` | CSV | Combined feature weights for all bootstrap iterations across 5 EUR:AFR composition ratios. Includes columns: `CpG`, `Coefficient`, `Ratio`, `Iteration`. |
| `Composition_summary.csv` | CSV | Performance metrics (R², MAE, RMSE, cross-validation metrics) across all bootstrap iterations and ratio gradients. |
| `full_clock_coef.csv` | CSV | Non-zero feature coefficients selected by the Elastic Net model trained on the full combined dataset. |
| `full_clock_summary.csv` | CSV | Final cross-validation fit and performance statistics for the full combined clock. |

---

### 3. Feature Characterization (`results/characterization/`)

**Process**: `CHARACTERIZE_CPGS`

Contains downstream characterization files assessing CpG stability, genomic locus overlap, delta-beta differences, and gene annotations across composition ratios.

| File Pattern | Format | Description |
|---|---|---|
| `cpg_characterization_*` | CSV / RDS | Detailed CpG-level characterization tables including annotated CpG lists, delta-beta stats, and genomic feature distributions across ratio models. |
| `core_genes_custom.csv` | CSV | Prioritized core gene list mapped from stable CpG markers consistently selected across bootstrap iterations. |

---

### 4. Clock Evaluation (`results/evaluation/`)

**Process**: `EVALUATE_CLOCKS`

Contains validation outputs from evaluating trained clocks on the independent MAGENTA validation cohort (GSE338167).

| File | Format | Description |
|---|---|---|
| `magenta_all_predicted_ages_and_acceleration.csv` | CSV | Sample-level predictions on the MAGENTA cohort. Columns include `Sample_ID`, `Chronological_Age`, `Predicted_Age`, `Age_Acceleration` (residuals), and `Ancestry`. |
| `composition_evaluation_summary.csv` | CSV | Aggregate performance metrics (Pearson correlation $r$, Spearman $\rho$, Median Absolute Error MAE, RMSE) evaluated per ancestry subgroup. |

---

### 5. Publication Figures (`results/figures/`)

**Process**: `GENERATE_FIG5`

Contains composite Figure 5 visualizations summarizing CpG feature stability, weight distributions, and clock performance across ancestry ratios.

| File | Format | Description |
|---|---|---|
| `Figure5.pdf` | Vector PDF | Publication-ready multi-panel vector graphic of Figure 5. |
| `Figure5.png` | High-Res PNG | High-resolution raster export of Figure 5 (300 DPI). |

---

### 6. Aim 3 AD Multivariable Integration Models (`results/aim3_models/`)

**Processes**: `LOGISTIC_MODELS`, `STRATIFIED_MODELS`, `ELASTIC_NET_METALEARNER`, `QUADRANT_ANALYSIS`

Contains statistical modeling results integrating predicted epigenetic age acceleration with Polygenic Risk Scores (PRS) and clinical covariates to assess Alzheimer's Disease (AD) risk.

| File | Format | Description |
|---|---|---|
| `primary_logistic_models.csv` | CSV | Multivariable logistic regression results for AD diagnosis, reporting Odds Ratios (OR), 95% Confidence Intervals, and p-values for age acceleration and covariates. |
| `stratified_logistic_models.csv` | CSV | Ancestry-stratified AD risk logistic model metrics (EUR-specific and AFR-specific associations). |
| `metalearner_feature_weights.csv` | CSV | Feature importance weights from the Elastic Net meta-learner integrating epigenetic clock predictions, PRS, and clinical risk factors. |
| `quadrant_concordance_stats.csv` | CSV | Concordance and risk enrichment statistics across the 4 risk quadrants defined by high/low epigenetic age acceleration and high/low genetic PRS risk. |
