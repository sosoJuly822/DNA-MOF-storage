# DNA-MOF-storage

This repository hosts research code and example processed results for a DNA storage study based on Metal-Organic Frameworks (MOFs).

## Project Overview

This project focuses on high-density DNA-based data storage using MOF-assisted preservation and sequencing-based data recovery. The repository provides data processing scripts for two major parts of the study:

- **Decoding workflow (clustering-based sequence generation)**: generating representative sequences for decoding from sequencing reads under the assumption that primer sequences, indices, and sequence segment lengths are known. The full decoding algorithm for reconstructing original stored information is described in a separate study (http://arxiv.org/abs/2601.16518) and is not included in this repository.
- **Error rate analysis workflow**: evaluating sequencing and decoding fidelity by calculating substitutions, insertions, and deletions.

The code is intended to support the data processing and analysis procedures used in our DNA-MOF storage study.

## Repository Structure

```text
DNA-MOF-storage/
├── README.md
├── scripts/
│   ├── decoding/                  # Scripts for clustering-based sequence preparation
│   │   ├── *.py                   # Preprocessing, index assignment, clustering, etc.
│   │   ├── *.ipynb                # Visualization and exploratory notebooks
│   │   └── *.sh                   # Batch processing / pipeline scripts
│   └── error_rate_analysis/       # Scripts for calculating sequencing and decoding errors
│       ├── *.py
│       ├── *.ipynb
│       └── *.sh
├── results/
│   ├── decoding/                  # Example processed sequences from clustering
│   └── error_rate_analysis/       # Example processed error rate analysis results
├── examples/                      # Small input/output examples for testing
└── docs/                          # Documentation of workflows (optional)
```

## Decoding Workflow (Clustering-based Sequence Preparation)

The decoding scripts located in `scripts/decoding/` are used to generate sequences suitable for downstream decoding. This workflow assumes:

- Known primer sequences
- Known index sequences
- Known length of each sequence segment

The clustering step organizes reads according to these constraints to produce representative sequences. **Full decoding of the stored digital information is not included in this repository.**

Example processed results are provided in `results/decoding/`.

## Error Rate Analysis Workflow

The error rate analysis scripts are located in `scripts/error_rate_analysis/`. This workflow evaluates sequencing and decoding accuracy:

- Matching sequencing reads or clustered sequences to reference sequences
- Identifying covered and uncovered sequences
- Counting substitutions, insertions, and deletions
- Calculating per-read, per-item, and overall error rates
- Summarizing error distributions

Example processed results are provided in `results/error_rate_analysis/`.

## Data Availability

Due to the large size of raw sequencing files and intermediate processing files, the complete raw datasets are not included. This repository currently provides core analysis scripts and representative processed results.

Full raw datasets will be made available through an appropriate public repository after publication. Before public release, data may be available from the corresponding authors upon reasonable request.

## Contact

For questions about code, data processing workflow, or access to raw data, please contact the corresponding authors or open an issue in this repository.

## Status

This repository is under active development. Code, example results, documentation, and data availability information will be updated during the manuscript submission and revision process.