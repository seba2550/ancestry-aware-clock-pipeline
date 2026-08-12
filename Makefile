# Makefile for Ancestry-Aware Methylation Clock Pipeline

.PHONY: help build-containers test clean

help:
	@echo "Available commands:"
	@echo "  make build-containers   - Build all Docker containers locally"
	@echo "  make test               - Run Nextflow pipeline with test profile"
	@echo "  make clean              - Clean Nextflow execution artifacts (.nextflow, work/)"

build-containers:
	@echo "Building Docker containers..."
	docker build -t sgonzalez/clock-training:1.0.0 containers/clock-training/
	docker build -t sgonzalez/r-bioconductor:1.0.0 containers/r-bioconductor/
	docker build -t sgonzalez/prs-tools:1.0.0 containers/prs-tools/

test:
	@echo "Running test profile..."
	nextflow run main.nf -profile test

clean:
	@echo "Cleaning execution logs and work directories..."
	rm -rf .nextflow* work/ results/
