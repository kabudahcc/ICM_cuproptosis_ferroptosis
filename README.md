# Cuproptosis- and Ferroptosis-Related Genes as Novel Biomarkers for Ischemic Cardiomyopathy

Analysis code for the manuscript:

> **Cuproptosis- and Ferroptosis-Related Genes as Novel Biomarkers for Ischemic Cardiomyopathy: An Integrated Bioinformatics and Machine Learning Study**
> (submitted to *Scientific Reports*)

## Overview

This repository contains the complete R analysis pipeline for an integrated
bioinformatics study of cuproptosis- and ferroptosis-related genes in ischemic
cardiomyopathy (ICM):

- Differential expression analysis of three discovery microarray cohorts
- Cross-dataset meta-analysis (ComBat + limma; permutation-based rank product, 8 datasets)
- GO/KEGG enrichment and compareCluster analysis
- Protein–protein interaction (PPI) network analysis (STRING)
- Weighted gene co-expression network analysis (WGCNA)
- Immune infiltration analysis (GSVA + MCPcounter) and immune checkpoint analysis
- Machine learning diagnostic models (LASSO-logistic regression, Random Forest, SVM)
  with 10-fold CV and leave-one-dataset-out cross-validation
- Unified external validation across six independent datasets
- Two-sample Mendelian randomization (eQTLGen cis-eQTLs → HERMES heart-failure GWAS)
- Single-cell RNA-seq validation (Seurat v5, GSE145154)

## Data sources

All datasets are publicly available:

| Dataset | Type | Use |
|---|---|---|
| GSE16499, GSE5406, GSE57338 | Microarray | Discovery cohorts |
| GSE116250, GSE55296, GSE42955, GSE26887, GSE46224, GSE52601 | RNA-seq / microarray | External validation |
| GSE145154 | scRNA-seq (10x) | Single-cell validation |
| FerrDb V3 (http://www.zhounan.org/ferrdb/) | Gene set | Ferroptosis genes |
| eQTLGen (2019-12-11 release) | cis-eQTL summary stats | MR exposure |
| HERMES GWAS (GWAS Catalog GCST009541) | GWAS summary stats | MR outcome |

## Pipeline (scripts/)

Run order approximately follows the step numbering:

| Script | Content |
|---|---|
| `step1_GEO_GSE16499.R`, `step1_GEO_GSE5406.R`, `step1_GEO_GSE57338.R` | GEO download, preprocessing, limma DEG analysis |
| `step2_venn_three.R`, `step2b_cell_death_DEGs.R`, `step2c_cell_death_heatmap.R` | Cross-dataset intersection, cuproptosis/ferroptosis DEGs |
| `step3b_cell_death_GO_KEGG.R`, `step4b_compareCluster_KEGG.R` | Functional enrichment |
| `step5_PPI_network.R`, `step5c_PPI_double_circle.R` | PPI network and hub genes |
| `step6_immune_infiltration.R`, `step11b_checkpoint_ICM.R`, `step11c_immune_cor_ICM.R` | Immune infiltration, checkpoints, correlations |
| `step7_machine_learning.R`, `step7b_LOOCV.R` | Diagnostic models (10-fold CV, LOOCV) |
| `step8_summary_and_figures.R` | Summary figures |
| `step9_external_validation.R`, `step9b_external_validation_GSE55296_GSE42955.R`, `step9c_external_validation_unified.R` | External validation |
| `step17_unified_validation_6datasets.R` | **Unified external validation, all six datasets (final version)** |
| `step10_single_cell.R` | Single-cell analysis (Seurat v5) |
| `step11a_WGCNA_ICM.R` – `step11a4_WGCNA_6panel_vector_ICM.R` | WGCNA and six-panel figure |
| `step12_pathway_scoring.R` | GSVA pathway scoring |
| `step13b_prepare_GSE52601.R`, `step13c_inspect_GSE141910.R` | Dataset preparation/inspection |
| `step14_meta_analysis_DEG.R` | Cross-dataset meta-analysis (ComBat + limma, rank product) |
| `step15_MR_analysis.R` | Two-sample Mendelian randomization |
| `step16_PCA_batch_effect.R` | PCA before/after ComBat (Supplementary Figure S6) |
| `extract_ferrdb.R` | FerrDb gene set extraction |

## Requirements

- R ≥ 4.4 (developed under R 4.6.0)
- CRAN: `ggplot2`, `dplyr`, `tidyr`, `glmnet`, `randomForest`, `e1071`, `pROC`,
  `caret`, `PRROC`, `igraph`, `pheatmap`, `FactoMineR`, `factoextra`, `ggrepel`, `data.table`
- Bioconductor: `GEOquery`, `limma`, `clusterProfiler`, `GSVA`, `sva`, `WGCNA`,
  `Seurat` (≥ 5.0), `edgeR`, `hugene10sttranscriptcluster.db`, `AnnotationDbi`

Raw data files (GEO series matrices, count tables) are **not** included and must
be downloaded from GEO before running the pipeline.

## License

MIT (see `LICENSE`).
