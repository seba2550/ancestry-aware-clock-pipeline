# Pipeline Parameters Reference

This document provides a comprehensive list of all configuration parameters for the **Ancestry-Aware Methylation Clock Nextflow Pipeline**, including default values, data types, and detailed functional descriptions.

Default parameter values are declared in [`params.yml`](file:///Users/sgonzalez/Desktop/Capra%20Lab/Thesis%20Project/ancestry-aware-clock-pipeline/params.yml) and [`nextflow.config`](file:///Users/sgonzalez/Desktop/Capra%20Lab/Thesis%20Project/ancestry-aware-clock-pipeline/nextflow.config).

---

## 1. Input Data & GEO Accessions

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--geo_eur_accession` | String | `"GSE55763"` | GEO accession ID for European training cohort dataset. |
| `--geo_afr_accession` | String | `"GSE210255"` | GEO accession ID for African training cohort dataset. |
| `--geo_magenta_accession` | String | `"GSE338167"` | GEO accession ID for MAGENTA validation cohort dataset. |
| `--eur_beta` | Path | `null` | Path to pre-processed European (EUR) beta-value matrix file. |
| `--eur_meta` | Path | `null` | Path to European (EUR) sample metadata file. |
| `--afr_beta` | Path | `null` | Path to pre-processed African (AFR) beta-value matrix file. |
| `--afr_meta` | Path | `null` | Path to African (AFR) sample metadata file. |
| `--magenta_beta` | Path | `null` | Path to MAGENTA validation cohort beta-value matrix file. |
| `--magenta_meta` | Path | `null` | Path to MAGENTA validation cohort sample metadata file. |

---

## 2. Output & Workflow Directives

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--outdir` | Path | `"./results"` | Target directory where all published pipeline outputs and figures are saved. |

---

## 3. Sample Filtering & Data Splitting

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--age_min` | Float | `19.0` | Minimum sample age threshold for inclusion in training/evaluation. |
| `--age_max` | Float | `78.0` | Maximum sample age threshold for inclusion in training/evaluation. |
| `--train_split_ratio` | Float | `0.80` | Proportion of data reserved for model training (vs. internal validation). |

---

## 4. Model Hyperparameters & Training Settings

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--alpha` | Float | `0.5` | Elastic Net mixing parameter $\alpha \in [0, 1]$. An $\alpha = 0.5$ balances L1 (Lasso) and L2 (Ridge) regularization penalties equally. |
| `--n_bootstrap` | Integer | `50` | Number of bootstrap iterations executed for each composition ratio gradient (100:0, 75:25, 50:50, 25:75, 0:100). |
| `--fixed_lambda` | Float | `0.010474566` | Pre-optimized regularization parameter $\lambda$ used for Elastic Net model fitting. |
| `--cv_folds` | Integer | `10` | Number of cross-validation folds used during internal model tuning. |

---

## 5. meQTL & Reference Annotation Assets

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--meqtl_cosmo` | Path | `"${projectDir}/assets/meqtl_lists/cosmo_meqtls.csv"` | Path to COSMO meQTL list for genomic characterization. |
| `--meqtl_genoa` | Path | `"${projectDir}/assets/meqtl_lists/genoa_meqtls.csv"` | Path to GENOA meQTL list for genomic characterization. |
| `--meqtl_epigen` | Path | `"${projectDir}/assets/meqtl_lists/epigen_meqtls.csv"` | Path to Epigen meQTL list for genomic characterization. |

---

## 6. Polygenic Risk Score (PRS) Weights

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--prs_weights` | Path | `"${projectDir}/assets/prs_weights/PGS003958_hmPOS_GRCh37.txt"` | Path to PGS Catalog weight file (PGS003958) used in Aim 3 multivariable integration models. |

---

## 7. Resource Limits & Execution Settings

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--max_memory` | Memory | `'32.GB'` | Maximum memory ceiling per task attempt. |
| `--max_cpus` | Integer | `8` | Maximum CPU thread ceiling per task attempt. |
| `--max_time` | Duration | `'24.h'` | Maximum time execution limit per task attempt. |
| `--container_registry` | String | `'docker.io/sgonzalez'` | Docker Hub registry prefix for pipeline container images. |
| `--slurm_queue` | String | `null` (defaults to `'standard'`) | SLURM partition / queue name when using `-profile slurm`. |

---

## 8. Test Mode Directives

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--is_test_run` | Boolean | `false` (set `true` in `-profile test`) | Flag indicating test run execution with synthetic data. |
| `--test_data_dir` | Path | `"${projectDir}/tests/test_data"` | Path to synthetic input datasets for test runs. |
