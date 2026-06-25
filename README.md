# Wheat GEA Framework

A reproducible pipeline for genome–environment association (GEA) analysis and adaptive germplasm prioritization in wheat.

## Overview

This project implements a modular workflow for identifying genomic variants associated with climatic variables while accounting for population structure and relatedness. Multiple genome-environment association methods and structure correction strategies are evaluated to identify robust climate-adaptive signals and prioritize germplasm accessions for downstream breeding and conservation applications.

## Pipeline Overview

1.  Environmental data extraction and processing

2.  Genotype quality control

3.  Population structure analysis

4.  Kinship analysis

5.  Genome-environment association analysis

    - LMM (GEMMA)
    - LFMM
    - RDA
    - pcadapt

6.  Evaluation and comparison of structure correction strategies

7.  Consensus SNP set construction

8.  Robustness-based consensus selection

9.  Accession-level adaptive germplasm scoring

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

These input file paths should be updated in the pipeline configuration file. 

## Configuration

Pipeline settings are controlled through:

``` text
config.yaml
```

This file contains:

- Input and output paths
- Command line tool paths
- Climate extraction settings
- Quality control thresholds
- Population structure parameters
- LMM, LFMM, RDA, and pcadapt settings
- Consensus analysis options

All configuration settings can be modified if needed; however, settings required by the current pipeline structure are annotated with `# DNC` (“Do Not Change”). In most cases, users should only modify settings that are not marked with `# DNC`, unless they are also updating the corresponding project file structure, column names, or pipeline code.


## Current Status

This repository is under active development.
