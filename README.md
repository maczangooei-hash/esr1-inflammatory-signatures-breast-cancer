# Inflammatory Signatures and Estrogen-Response Programs in Breast Cancer

## Study

**Inflammatory Transcriptional Signatures Are Associated with Estrogen-Response Programs Beyond ESR1 mRNA Abundance in ER-Positive Breast Cancer**

This repository contains selected analysis scripts, results, and Figure 1 for the study of associations between inflammatory Hallmark scores and estrogen-response Hallmark scores in ER-positive breast tumors.

## Data sources

- GSE96058: clinical breast-cancer expression data.
- GSE126003: exploratory T47D siRNA experiment.

The two platform-defined subsets of GSE96058 belong to the same GEO series and are not independent external validation cohorts.

Download source data from the [NCBI Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/). Raw datasets, large expression matrices, and the Hallmark gene-set file are not included.

## Contents

- `scripts/clinical/`: clinical analysis and sensitivity scripts.
- `scripts/experimental/`: exploratory GSE126003 analysis.
- `results/manuscript_tables/`: selected manuscript tables.
- `results/clinical/`: selected supporting results.
- `figures/`: Figure 1 and its generation script.

## Software

The R analyses were developed with R 4.5.0. Scripts using `limma` require the Bioconductor package `limma`.

Figure generation requires Python, NumPy, and Matplotlib. From the repository root, run:

    python figures/create_figure1_primary.py

The script generates PDF and 600-dpi TIFF versions of Figure 1.

## Reproducibility

The complete R workflow is not yet reproducible from a fresh clone alone because required source data and intermediate files are not included. Input paths and sample/gene alignment must be checked before running the analyses.

Hallmark scores represent gene-expression patterns, not direct measurements of estrogen-receptor protein activity. The analyses establish associations, not causation. The GSE126003 analysis uses processed RPKM transformed as `log2(RPKM + 1)`, not raw counts.

## Data availability

Source datasets are publicly available from GEO under accessions GSE96058 and GSE126003. The public repository URL will be added after release.

## License

A license has not yet been selected. Redistribution terms for included materials must be checked before public release.
