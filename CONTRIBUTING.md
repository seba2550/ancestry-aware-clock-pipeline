# Contributing to Ancestry-Aware Epigenetic Clock Pipeline

Thank you for your interest in contributing to the **Ancestry-Aware Epigenetic Clock Pipeline**! We welcome contributions from bioinformatics researchers, computational biologists, software engineers, and epigenetics developers.

This document outlines the guidelines and standards for contributing to this repository.

---

## Table of Contents

1. [Code of Conduct](#code-of-conduct)
2. [Getting Started](#getting-started)
3. [Development Environment](#development-environment)
4. [Pipeline Architecture & Conventions](#pipeline-architecture--conventions)
   - [Nextflow Standards (DSL2)](#nextflow-standards-dsl2)
   - [Python Guidelines](#python-guidelines)
   - [R Guidelines](#r-guidelines)
   - [Containerization Standards](#containerization-standards)
5. [Testing & Quality Assurance](#testing--quality-assurance)
6. [Submitting Contributions](#submitting-contributions)
   - [Issue Reporting](#issue-reporting)
   - [Commit Conventions](#commit-conventions)
   - [Pull Request Checklist](#pull-request-checklist)

---

## Code of Conduct

We are committed to fostering an inclusive, respectful, and welcoming scientific community. Please ensure all interactions—whether in issues, pull requests, or discussions—remain professional and constructive.

---

## Getting Started

1. **Fork & Clone** the repository:
   ```bash
   git clone https://github.com/YOUR_USERNAME/ancestry-aware-clock-pipeline.git
   cd ancestry-aware-clock-pipeline
   ```
2. **Set up a feature branch**:
   ```bash
   git checkout -b feature/my-new-feature
   ```

---

## Development Environment

To develop and test pipeline modules locally, ensure you have the following installed:

- **Nextflow**: `>=22.10.0`
- **Container Engine**: Docker (`>=20.10`) or Singularity (`>=3.8`)
- **Python**: `>=3.10` (for local script execution and testing)
- **R / Bioconductor**: `>=4.2` (for downstream statistical analysis and plotting)
- **Make**: Standard GNU Make utility

Verify your environment using:
```bash
nextflow info
docker --version
```

---

## Pipeline Architecture & Conventions

### Nextflow Standards (DSL2)

- **Modularity**: All process logic must reside in discrete DSL2 module files under `modules/local/` or `modules/nf-core/`.
- **Process Labels**: Always specify resource and container labels (e.g., `label 'process_medium'`, `label 'container_python'`). Avoid hardcoding CPU/memory limits directly within process definitions.
- **Parameters**: Define default parameter values in `params.yml` or `nextflow.config`. Override defaults using profile configs in `conf/`.
- **Clean Inputs/Outputs**: Use explicit file patterns and emitted channels with named outputs for pipeline composability.

### Python Guidelines

- **Code Style**: Adhere to PEP 8 syntax formatting. Use `black` and `flake8` for formatting and linting.
- **Type Annotations**: Provide type hints for public functions and classes wherever practical.
- **CLI Interface**: Python scripts in `bin/` must use `argparse` or `click` to handle command-line options cleanly and produce informative `--help` outputs.
- **Dependencies**: Keep dependencies lightweight and pin exact version specs in environment YAML files (`containers/environment.yml`).

### R Guidelines

- **Style**: Follow tidyverse style guidelines (consistent indentation, clear variable names).
- **Bioconductor**: Use standard Bioconductor structures (`SummarizedExperiment`, `GRanges`) for epigenomics data manipulation.
- **Error Handling**: Use explicit error checking and structured exits (`stop()` with diagnostic messages) in `bin/` scripts.

### Containerization Standards

- Each pipeline process label maps to a pre-built container image specified in `nextflow.config`.
- Dockerfiles located in `containers/` must:
  - Be built from minimal base images (e.g., `python:3.10-slim`, `rocker/r-ver:4.2.2`).
  - Pin dependency versions explicitly.
  - Set appropriate file permissions and execution non-root user settings.

---

## Testing & Quality Assurance

Before submitting changes, run the test suite to verify pipeline execution and code formatting:

1. **Run Integration Test Profile**:
   ```bash
   nextflow run main.nf -profile test,docker
   ```
2. **Execute Makefile Verification**:
   ```bash
   make test
   make lint
   ```
3. **Validate Configuration**:
   Ensure configuration syntax is clean and all input profile paths resolve properly.

---

## Submitting Contributions

### Issue Reporting

- **Bug Reports**: Open an issue describing the expected vs. actual behavior, complete with Nextflow execution logs (`.nextflow.log`) and environment info.
- **Feature Requests**: Describe the biological rationale, algorithm, or pipeline stage you wish to add.

### Commit Conventions

We follow Conventional Commits for clear versioning history:
- `feat:` A new pipeline feature or module
- `fix:` A bug fix in pipeline script or module execution
- `docs:` Documentation updates
- `refactor:` Code improvements that do not alter output logic
- `test:` Adding or modifying integration tests or sample test datasets
- `ci:` Updates to GitHub Actions workflows

Example:
```bash
git commit -m "feat(modules): add ancestry-stratified ElasticNet training module"
```

### Pull Request Checklist

When opening a Pull Request (PR):

- [ ] My code follows the code style guidelines of this project.
- [ ] I have executed `nextflow run main.nf -profile test,docker` and verified that all test tasks succeed.
- [ ] I have updated relevant documentation in `docs/` and inline comments.
- [ ] I have added appropriate process labels and container specs for any new processes.
- [ ] My commit messages adhere to Conventional Commits format.

Thank you for contributing to open and reproducible epigenetics research!
