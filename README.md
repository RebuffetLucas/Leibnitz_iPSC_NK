# Analysis of NK cells products derived from a new cell product manufacturing protocol

## Project Overview

This project focuses on understanding the transcriptional and functional diversity of K cells products derived from a new cell product manufacturing protocol. It employs state-of-the-art single-cell RNA sequencing (scRNA-seq) and computational methods to:

1. Identify the composition of the populations derived from the NK progenitors.
2. Map them to some of the current NK cell classification, using methods such as Transfer mapping or evaluation of the overlap in the signatures
3. Assess the functions of the NK cells produced by this system


---

## Key Analyses

### 1. **Data Preprocessing and Integration**
- **Preprocessing:**
  - Normalized and scaled datasets for FL, BM, and small intestine (SI).
  - Identified variable features using `Seurat`
- **Automatic identification of NK cells and other immune cells:**
  - With SingleR
  - With Azimuth ( a reference based approach).

### 2. **NK Subset Mapping**
- We try to map the NK1, NK2, NKint and NK3 populations using transfer mapping approach:
- We verify subset-specific signatures using DotPlots and FeaturePlots.


### 3. **Visualization**
- PCA and UMAP for dimensionality reduction.
- DotPlots and FeaturePlots for gene expression.
- Heatmaps for regulon activity and pseudotime-ordered gene expression.
- Barplots for proportional comparisons across tissues and clusters.

---

## How to Reproduce the Analysis

1. **Install Dependencies:**
   - R packages: `Seurat`, `Harmony`, `SCENIC`, `ggplot2`, `ComplexHeatmap`, etc.
   - Ensure all required scripts and raw data are available in the repository.

2. **Run Scripts:**
   - Follow the order in the `03_Scripts/` directory to reproduce each step of the analysis.

3. **Generate Figures:**
   - Visualization scripts automatically save outputs to the `Figures/` directory.

---

## Acknowledgments

For questions or issues, contact [rebuffet@ciml.univ-mrs.fr]


