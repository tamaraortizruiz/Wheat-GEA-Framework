# Wheat GEA Framework

A reproducible pipeline for genome–environment association (GEA) analysis and adaptive germplasm prioritization in wheat.

## Overview

This project implements a modular workflow for identifying genomic variants associated with climatic variables while accounting for population structure and relatedness. Multiple genome-environment association methods and structure correction strategies are evaluated to identify robust climate-adaptive signals and prioritize germplasm accessions for downstream breeding and conservation applications.

## Pipeline Overview

1.  Quality control of genotype data using PLINK

2.  Climate variable extraction from WorldClim

3.  Environmental filtering

4.  Population structure analysis using PCA

5.  Kinship matrix calculation using GEMMA

6.  Genome–environment association analysis

    - LMM (GEMMA)
    - LFMM
    - RDA
    - pcadapt

7.  Evaluation and comparison of structure correction strategies

8.  Consensus SNP identification

9.  Robustness-based consensus selection

10. Accession-level adaptive germplasm scoring

## Project Structure

``` text
RawData/     Input genotype and metadata files (not tracked by Git)
Data/        Intermediate files generated during analysis (not tracked by Git)
Output/      Analysis results, figures, and summary tables (not tracked by Git)
Functions/   Modular R functions used throughout the pipeline
```

## Required Input Files

The following files should be placed in `RawData/` before running the pipeline:

- PLINK genotype files (`.bed`, `.bim`, `.fam`)
- Passport metadata containing sample identifiers and geographic coordinates (longitude and latitude)
- Optional sample subsets or marker subsets

## Configuration

Pipeline settings are controlled through:

``` text
config.yaml
```

This file contains:

- Input and output paths
- Command line tool paths
- Quality control thresholds
- Climate extraction settings
- Population structure parameters
- LMM, LFMM, RDA, and pcadapt settings
- Consensus analysis options

## Current Status

This repository is under active development.
