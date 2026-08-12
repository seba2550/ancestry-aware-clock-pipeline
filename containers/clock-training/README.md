# Clock Training Container (Python 3.11)

This repository contains the container configuration for the Python-based epigenetic clock training environment in the ancestry-aware clock pipeline.

## Directory
`/Users/sgonzalez/Desktop/Capra Lab/Thesis Project/ancestry-aware-clock-pipeline/containers/clock-training/`

## Container Purpose
This container packages Python 3.11 with core scientific, statistical, and machine learning libraries required for penalised regression modeling (ElasticNet, GLMNet) and training ancestry-aware epigenetic clock models.

## Dependencies

### System Dependencies
- `build-essential` (C/C++ compiler toolchain)
- `gfortran` (Fortran compiler required for compiling numerical routines like `glmnet`)
- `libopenblas-dev` (Optimised BLAS routines)
- `git` (Version control integration)

### Python Libraries
- `numpy`
- `pandas`
- `scikit-learn`
- `scipy`
- `glmnet`
- `joblib`

## Build Command

To build the Docker image with tag `sgonzalez/clock-training:1.0.0`:

```bash
docker build -t sgonzalez/clock-training:1.0.0 .
```

## Validation & Execution Instructions

To execute the validation script within the container:

```bash
docker run --rm sgonzalez/clock-training:1.0.0 python test_container.py
```

To launch an interactive session inside the container:

```bash
docker run --rm -it -v $(pwd):/app sgonzalez/clock-training:1.0.0 /bin/bash
```

## Included Files
- `Dockerfile`: Multi-stage / minimal Debian-slim based container definition.
- `test_container.py`: Standalone Python validation script verifying dependency imports and running a 10-sample ElasticNet regression.
- `README.md`: Container documentation and instructions.
