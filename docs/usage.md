# Pipeline Usage Guide

This document provides detailed execution instructions for the **Ancestry-Aware Methylation Clock Nextflow Pipeline**, covering local container execution (Docker/Singularity), HPC SLURM cluster execution, test mode execution, and parameter customization.

---

## Prerequisites

Before running the pipeline, ensure the following software is installed on your system:

- **Nextflow**: Version `!>=22.10.0`
- **Java**: Java 11 or later (required by Nextflow)
- **Container Engine** (one of the following):
  - **Docker**: Desktop or Engine (`>=20.10`) for local runs.
  - **Singularity / Apptainer**: (`>=3.8`) for HPC environments.

---

## Execution Modes

### 1. Test Mode Execution

Test mode uses minimal synthetic datasets packaged within `tests/test_data/` to validate pipeline syntax, container execution, and process orchestration. It automatically sets `--n_bootstrap 2` and reduces resource limits.

#### Local Execution (Synthetic Data Only)
```bash
nextflow run main.nf -profile test
```

#### Local Execution with Docker
```bash
nextflow run main.nf -profile test,docker
```

#### Local Execution with Singularity
```bash
nextflow run main.nf -profile test,singularity
```

#### Via Makefile
```bash
make test
```

---

### 2. Local Production Execution (Docker)

For full pipeline execution on a local machine or workstation using Docker:

```bash
nextflow run main.nf \
  -profile docker \
  --eur_beta /path/to/GSE55763_beta.rds \
  --eur_meta /path/to/GSE55763_meta.rds \
  --afr_beta /path/to/GSE210255_beta.rds \
  --afr_meta /path/to/GSE210255_meta.rds \
  --magenta_beta /path/to/GSE338167_beta.rds \
  --magenta_meta /path/to/GSE338167_meta.rds \
  --outdir ./results
```

---

### 3. Local Production Execution (Singularity)

If executing locally using Singularity/Apptainer instead of Docker:

```bash
nextflow run main.nf \
  -profile singularity \
  --eur_beta /path/to/GSE55763_beta.rds \
  --eur_meta /path/to/GSE55763_meta.rds \
  --afr_beta /path/to/GSE210255_beta.rds \
  --afr_meta /path/to/GSE210255_meta.rds \
  --magenta_beta /path/to/GSE338167_beta.rds \
  --magenta_meta /path/to/GSE338167_meta.rds \
  --outdir ./results
```

---

### 4. HPC Cluster Execution (SLURM + Singularity)

For high-performance computing clusters managed by the SLURM workload manager:

```bash
nextflow run main.nf \
  -profile slurm,singularity \
  --slurm_queue standard \
  --eur_beta /scratch/data/GSE55763_beta.rds \
  --eur_meta /scratch/data/GSE55763_meta.rds \
  --afr_beta /scratch/data/GSE210255_beta.rds \
  --afr_meta /scratch/data/GSE210255_meta.rds \
  --magenta_beta /scratch/data/GSE338167_beta.rds \
  --magenta_meta /scratch/data/GSE338167_meta.rds \
  --outdir /scratch/user/clock_results
```

#### Customizing Resource Constraints for Cluster Runs
You can adjust maximum host hardware boundaries directly from the command line:

```bash
nextflow run main.nf \
  -profile slurm,singularity \
  --slurm_queue gpu_long \
  --max_memory '64.GB' \
  --max_cpus 16 \
  --max_time '48.h' \
  --outdir /scratch/user/clock_results
```

---

## Parameter Override Examples

Parameters can be overridden via command-line arguments (`--parameter_name value`) or custom parameter configuration files (`-params-file my_params.yml`).

### Command-Line Overrides

#### Adjusting Bootstrap Iterations and Elastic Net Alpha
```bash
nextflow run main.nf \
  -profile docker \
  --n_bootstrap 100 \
  --alpha 0.75 \
  --fixed_lambda 0.015 \
  --outdir ./results_custom_alpha
```

#### Downsampling Age Range Filter
```bash
nextflow run main.nf \
  -profile docker \
  --age_min 30.0 \
  --age_max 65.0 \
  --train_split_ratio 0.85
```

### Using a Custom Parameters File

You can create a custom `custom_params.yml` file:

```yaml
n_bootstrap: 100
alpha: 0.5
fixed_lambda: 0.010474566337054866
outdir: "./results_production_run"
age_min: 20.0
age_max: 75.0
```

And execute Nextflow with:

```bash
nextflow run main.nf -profile docker -params-file custom_params.yml
```

---

## Execution Lifecycle & Workflow Options

### Resuming Failed or Interrupted Runs
Nextflow caches process execution state in `work/`. If a run is interrupted or fails due to transient node errors, resume it without re-executing finished tasks using `-resume`:

```bash
nextflow run main.nf -profile docker -resume
```

### Cleaning Up Temporary Files
To clear cached Nextflow execution artifacts (`.nextflow*`, `work/`, `results/`):

```bash
make clean
```
Or manually:
```bash
rm -rf .nextflow* work/ results/
```
