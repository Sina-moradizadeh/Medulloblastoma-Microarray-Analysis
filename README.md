# Medulloblastoma-Microarray-Analysis
A comprehensive microarray gene expression analysis pipeline for medulloblastoma using R, including quality control, annotation, differential expression analysis with limma, visualization, and functional enrichment analysis.
# Medulloblastoma Microarray Analysis

## Overview

This project presents a complete bioinformatics workflow for microarray gene expression analysis of medulloblastoma using the publicly available GEO dataset GSE124814.

The analysis was performed in R using Bioconductor packages and includes differential gene expression analysis, data visualization, and functional enrichment analysis.

This project was developed for educational purposes to practice reproducible bioinformatics analysis.

---

## Dataset

- GEO Accession: GSE124814
- Platform: Affymetrix Microarray
- Disease: Medulloblastoma

Note:
The expression matrix downloaded from GEO had already been preprocessed by the original study, including background correction, normalization, cross-platform integration, and batch-effect correction.

---

## Workflow

1. Download data from GEO
2. Data exploration
3. Quality control
4. Gene annotation
5. Differential expression analysis using limma
6. Volcano plot
7. Heatmap
8. GO enrichment analysis
9. KEGG pathway analysis

---

## Tools

- R
- GEOquery
- tidyverse
- limma
- AnnotationDbi
- org.Hs.eg.db
- clusterProfiler
- enrichplot
- ggplot2

---

## Repository Structure

├── 01_Download_GSE_Dataset.R
├── 02_Primary_QC.R
├── 03_Preprocessing.R
├── 04_DEG_Limma.R
├── 05_Visualization.R
├── 06_Enrichment_Analysis.R
├── Results/
└── README.md

---

## Results

The analysis identified differentially expressed genes between the selected biological groups. Functional enrichment analyses were performed to investigate significantly enriched biological processes and signaling pathways.

---

## Future Work

- Machine learning-based biomarker identification
- External dataset validation
- Network analysis
- Survival analysis

---

## Author

Sina Moradizadeh
