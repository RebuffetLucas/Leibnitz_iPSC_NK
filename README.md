# Single-Cell RNA-seq Analysis with Seurat and SingleR

This repository contains a modular and reproducible R-based workflow for processing and analyzing single-cell RNA sequencing (scRNA-seq) data using **Seurat** and **SingleR**. The pipeline performs quality control, dimensionality reduction, clustering, automated cell annotation, gene module scoring, and exploratory analyses such as metadata correlations.

## 🧩 Project Structure

- **`analysisParams.R`**: Centralized configuration file containing all parameter settings used throughout the analysis (e.g. thresholds, resolution, gene counts, PCA/UMAP dimensions, etc.).
- **`report_01_QC.Rmd`**: Main RMarkdown file that compiles the analysis report. This document sources all code chunks (labeled with `@knitr`) from external scripts for clean modular execution and documentation.

## 📁 Script Overview

### 1. `violin_QC`, `FeatureScatter`
Visualize core QC metrics:
- Violin plots of `nFeature_RNA`, `nCount_RNA`, `percent.mito`, and `percent.ribo` across samples
- Scatter plots to detect low-quality cells, doublets, and extreme outliers

### 2. `Variable_Features`, `scaling_and_pca`
Prepares data for dimensionality reduction and clustering:
- Log-normalization and selection of highly variable genes (`vst` method)
- Scaling and PCA computation
- Visual inspection of PCA embeddings

### 3. `UMAP`
Performs:
- Graph-based neighbor detection
- Clustering (with user-defined resolution)
- UMAP embedding for 2D visualization of clusters

### 4. `markers`
Detects cluster-specific markers:
- Finds positive markers using `FindAllMarkers`
- Filters by log fold-change and FDR
- Displays interactive tables with cluster-specific color coding

### 5. `FeaturePlot`
Visualizes gene expression:
- QC genes and NK-related markers
- Colored by expression intensity across the UMAP embedding

### 6. `Highlight_cells_to_remove`
Flags and highlights unwanted clusters:
- Visualizes suspect clusters (e.g. myeloid contamination) on UMAP
- Outputs barcodes of selected cells in a downloadable format

### 7. `Scoring_NK123_and_13Genes`
Computes gene module scores:
- NK1/NK2/NK3 signatures derived from external marker sets
- Crinier et al. (2018) 13-gene NK signature (loaded from Excel)
- Displays UMAP expression maps and violin plots of module scores

### 8. `SingleR`, `CellsToRemove`
Automated cell annotation:
- Uses **SingleR** with the **BlueprintEncodeData** reference for immune cell typing
- Stores predicted fine-level labels in the Seurat metadata
- Visualizes cell identity on UMAP, barplots cluster composition
- Highlights and exports barcodes for unwanted cell types (e.g. `MEP`, `GMP`, `Monocytes`)

### 9. `correlation_meta_data_PC`
Explores relationships between metadata and principal components:
- Computes Spearman correlation between PCs and numeric metadata variables
- Visualizes as a correlation heatmap using `corrplot`

## 🔧 Requirements

Main R packages:
- `Seurat`, `SingleR`, `celldex`
- `ggplot2`, `ggdist`, `gghalves`, `viridis`, `RColorBrewer`
- `openxlsx`, `readxl`
- `DT`, `htmltools`
- `dplyr`, `tidyr`, `corrplot`

## 📤 Outputs

- UMAP and violin plots
- Cluster marker tables (interactive and downloadable)
- Barcode lists for low-quality or unwanted cell populations
- Gene module score visualizations
- PCA-metadata correlation heatmap

## ⚙️ Parameters

All user-defined parameters — including thresholds for filtering, number of features, clustering resolution, and visualization settings — are defined in `analysisParams.R`. This ensures consistency and facilitates reproducibility across reruns.

## 📄 Report Generation

The full analysis pipeline is compiled and rendered from the RMarkdown file `report_01_QC.Rmd`. All code chunks are modular and sourced via knitr (`## @knitr`) to allow for clean organization and reusability.

## 👥 Authors

This pipeline was developed for immunological single-cell profiling projects, with a particular focus on NK cell heterogeneity and annotation.

---

Feel free to open an issue or PR to suggest improvements or extensions to the pipeline.