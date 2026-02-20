# Rscripts
# ER stress–autophagy intersectome in PCOS (analysis code & qPCR source data)

This repository contains the R scripts and source data used to reproduce the main analyses and figures in the manuscript:

**“An ER stress–autophagy intersectome identifies candidate diagnostic hub genes and potential microenvironmental interactions in PCOS”**

The study integrates public GEO datasets (bulk transcriptomics and scRNA-seq) with qPCR validation in human granulosa cells to identify and validate hub genes linked to ER stress and autophagy programs in PCOS.

---

## Contents

### 1) Code
- R scripts for:
  - GEO data download and preprocessing (probe annotation, log2 transform, normalization)
  - Batch correction (ComBat)
  - Differential expression analysis (limma)
  - Functional enrichment analysis (GO/KEGG)
  - Machine-learning feature selection and model development
  - Model evaluation (ROC, calibration, decision curve analysis)
  - Single-gene GSEA
  - Immune signature scoring (ssGSEA)
  - Figure generation

### 2) Data
- **qPCR source data** (Ct values) used for plotting and statistical testing.
- Public datasets are not redistributed here and should be obtained from GEO using the accession numbers listed below.

---

## Software environment

- **R version:** R software **v4.3.3**
- Recommended OS: Windows/macOS/Linux (tested under standard RStudio workflows)

### Key R packages (typical)
The scripts may use (depending on the module):
- GEOquery, limma, sva, ggplot2, pheatmap
- clusterProfiler, org.Hs.eg.db (or org.Mm.eg.db if needed), enrichplot
- GSVA, pROC, rms, rmda (or equivalent DCA package)
- caret / randomForest / e1071 / xgboost (for ML modules)
- Seurat (v5+) and SingleR (for scRNA-seq module)

> If you prefer fully reproducible dependency locking, consider running `renv::init()` after cloning the repo.

---

## Public datasets analysed (GEO)

Bulk granulosa cell expression datasets:
- **GSE34526**
- **GSE10946**
- **GSE106724**
- **GSE98595**

Single-cell dataset (mouse ovary):
- **GSE268919**

All bulk and scRNA-seq datasets are publicly available at NCBI GEO and can be downloaded using `GEOquery` or via the GEO website.
