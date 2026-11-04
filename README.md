# DNA-MOF-storage

This repository hosts research code and example processed results for a DNA storage study based on Metal-Organic Frameworks (MOFs).

The corresponding raw Illumina sequencing datasets are publicly available on  **Zenodo** :

[**https://zenodo.org/records/23938446**](https://zenodo.org/records/23938446)

## Project Overview

This project focuses on high-density DNA-based data storage using MOF-assisted preservation and sequencing-based data recovery. The repository provides data processing scripts for two major parts of the study:

- **Clustering and consensus workflow** : grouping sequencing reads according to known primer sequences, indices, and sequence segment lengths, followed by generation of representative consensus sequences for downstream analysis. The full decoding algorithm for reconstructing the original stored information is described in a separate study ([http://arxiv.org/abs/2601.16518](http://arxiv.org/abs/2601.16518)) and is not included in this repository.
- **Reference-based error rate analysis workflow** : matching sequencing reads to the known designed sequence library and quantifying sequencing errors, including substitutions, insertions, and deletions.

The repository also includes `refLib-SEQ2210.xlsx`, which contains the encoded DNA sequence library used in this study and serves as the reference sequence library for downstream sequencing data analysis.

The code is intended to support the data processing and analysis procedures used in our DNA-MOF storage study.

## Repository Structure

```text
DNA-MOF-storage/
├── README.md
├── refLib-SEQ2210.xlsx            # Encoded DNA sequence library used in this study
├── scripts/
│   ├── clustering_consensus/                  # Scripts for clustering-based sequence preparation
│   │   ├── *.py                   # Preprocessing, index assignment, clustering, etc.
│   │   ├── *.ipynb                # Visualization and exploratory notebooks
│   │   └── *.sh                   # Batch processing / pipeline scripts
│   ├── error_rate_analysis/       # Scripts for reference-based sequencing error analysis
│   │   ├── *.py
│   │   ├── *.ipynb
│   │   └── *.sh
|   └── downsampling/                  # Sequencing-depth downsampling
├── results/
│   ├── clustering_consensus/                  # Example processed sequences from clustering
│   └── error_rate_analysis/       # Example processed error rate analysis results
├── examples/                      # Small input/output examples for testing
└── docs/                          # Documentation of workflows (optional)
```

## Reference Sequence Library

The file `refLib-SEQ2210.xlsx`, located in the repository root at the same level as the `scripts/` and `results/` directories, contains the encoded DNA sequence library used in this study. These sequences represent the designed sequence library after encoding and are used as reference sequences in downstream sequencing data processing and error rate analysis.

## Clustering and Consensus Workflow

The scripts located in `scripts/clustering_consensus/` are used to classify sequencing reads and generate consensus sequences for downstream decoding. This workflow assumes:

* Known primer sequences
* Known index sequences
* Known length of each sequence segment

Reads are first assigned according to the known primer and index information. Within each assigned group, consensus sequences are generated from the sequencing reads without using the designed payload sequence as a reference.

The designed reference library may subsequently be used to evaluate the accuracy and recovery performance of the generated consensus sequences, but it is not used to construct the consensus payload sequence itself.

**The full decoding algorithm used to reconstruct the original stored digital information is not included in this repository.**

Example processed results are provided in `results/clustering_consensus/`.

## Reference-based Error Rate Analysis

The scripts located in `scripts/error_rate_analysis/` are used to characterize sequencing errors by direct comparison with the known designed sequence library.

Sequencing reads are assigned to their corresponding designed sequences using primer, index, and sequence-matching information. Each assigned read is then compared with its reference sequence to:

* Identify correctly sequenced positions
* Count substitutions, insertions, and deletions
* Calculate per-position and overall sequencing error rates
* Summarize error distributions across samples and sequences

The encoded sequence library provided in `refLib-SEQ2210.xlsx` serves as the reference library for this analysis.

Example processed results are provided in `results/error_rate_analysis/`.

## Data Availability

The raw Illumina sequencing datasets associated with this study are publicly available on  **Zenodo** :

[**https://zenodo.org/records/23938446**](https://zenodo.org/records/23938446)

The complete designed DNA sequence library used for downstream analysis is provided in this repository as `refLib-SEQ2210.xlsx`.

## Contact

For questions about the code or data processing workflow, please contact the corresponding authors or open an issue in this repository.

## Status

This repository accompanies the manuscript **“Interwoven DNA-MOF crystals for multi-millennial data storage.”** Code and documentation may be updated to improve reproducibility and usability.
