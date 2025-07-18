## @knitr FeaturePlot

# Print section header for RMarkdown/knitr
cat(" \n \n")
cat("#### FeaturePlot ")
cat(" \n \n")

# ---- FeaturePlot: QC metrics ----
# Visualize per-cell metadata (mito %, ribo %, gene/UMI counts) on UMAP
print(
  FeaturePlot(
    myObjectSeurat,
    features = c("percent.mito", "percent.ribo", "nFeature_RNA", "nCount_RNA"),
    pt.size = PT_SIZE_FEATURE_PLOT
  ) & 
    scale_color_viridis_c(option = "plasma")  # Use Viridis "plasma" colormap
)

# ---- FeaturePlot: Selected genes (first group) ----
# Visualize expression of immune-related or NK-associated genes
print(
  FeaturePlot(
    myObjectSeurat,
    features = LIST_GENE_1,
    pt.size = PT_SIZE_FEATURE_PLOT
  ) &
    scale_color_viridis_c(option = "plasma")
)

# ---- FeaturePlot: Alternative gene set (second group) ----
# Repeat similar plot with a different gene list
print(
  FeaturePlot(
    myObjectSeurat,
    features = LIST_GENE_2,
    pt.size = PT_SIZE_FEATURE_PLOT
  ) &
    scale_color_viridis_c(option = "plasma")
)





## @knitr Highlight_cells_to_remove

# ---- Highlight cells to remove ----
cat("#####" ,"UMAP","\n")

# Identify cells in cluster 3
cells <- WhichCells(myObjectSeurat, idents = 3)

# Plot UMAP with cluster 3 highlighted
print(
  DimPlot(
    myObjectSeurat,
    reduction = "umap",
    cells.highlight = cells,
    pt.size = UMAP_PT_SIZE,
    sizes.highlight = 1.5
  ) +
    scale_color_manual(
      labels = c("others", "cells to remove"),
      values = c("grey", "#CE4C4B")  # Red highlight for flagged cells
    )
)

cat(' \n \n')

# ---- Show barcodes of cells to remove ----
cat("#####" ,"Barcodes","\n")

# Display barcodes in a downloadable DataTable
DT::datatable(
  as.data.frame(cells),
  rownames = FALSE,
  extensions = 'Buttons',
  options = list(
    dom = 'Blfrtip',
    buttons = c('excel', "csv"),
    fixedHeader = TRUE
  )
)

cat(' \n \n')




## @knitr Scoring_NK123_and_13Genes

cat(" \n \n")
cat("#### Scoring NK 1, NK2, NK3 and 13 Genes Crinier ")
cat(" \n \n")

# ---- Load NK1/NK2/NK3 marker genes (from MetaNK dataset) ----
Markers_Seurat <- readRDS(file.path(PATH_EXPERIMENT_REFERENCE, "Tables_Gene_Signatures/All_Markers_MetaNK_CITEseq3clusters.rds"))

# ---- Filter top upregulated markers per cluster ----
top_All <- Markers_Seurat %>%
  filter(avg_log2FC > 0) %>%
  filter(p_val_adj < MARKERS_PVALUE_TRESHOLD) %>%
  filter(!grepl("RPL|RPS|MT-", gene)) %>%  # Remove ribosomal/mitochondrial genes for plotting
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_TOP_SCORING, wt = avg_log2FC)  # Select top markers by FC

# ---- Build list of gene sets per cluster ----
list_Top_Genes <- list()
for (i in levels(top_All$cluster)) {
  top_clust <- top_All %>% filter(cluster == i)
  list_Top_Genes <- append(list_Top_Genes, list(top_clust$gene))
}
names(list_Top_Genes) <- levels(top_All$cluster)

# ---- Format for AddModuleScore ----
List_To_Use <- lapply(list_Top_Genes, function(x) as.data.frame(x))
MONITORED_Markers <- List_To_Use

# ---- Score each NK module ----
# AddModuleScore creates a new metadata column for each gene set
for (i in names(MONITORED_Markers)) {
  myObjectSeurat <- AddModuleScore(
    myObjectSeurat,
    features = as.list(MONITORED_Markers[[i]]),
    pool = NULL,
    name = i,
    seed = SEED
  )
}

# ---- Plot scored modules ----
p2 <- DimPlot(myObjectSeurat, reduction = "umap")

# Display spatial expression on UMAP
p3 <- FeaturePlot(
  myObjectSeurat,
  features = paste0(names(MONITORED_Markers), "1"),  # "1" suffix added by AddModuleScore
  max.cutoff = MAX_CUTOFF_SCORE_MODULE
) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

# Display score distributions per cluster
p4 <- VlnPlot(
  myObjectSeurat,
  features = paste0(names(MONITORED_Markers), "1"),
  y.max = MAX_CUTOFF_SCORE_MODULE,
  pt.size = 0,
  ncol = 1
) & stat_summary(fun.data = data_summary, color = "black")

cat(' \n \n')
print(p3 + p2)  # Combine spatial and clustering plots
cat(' \n \n')

cat(' \n \n')
print(p4)       # Show violin plots of module scores
cat(' \n \n')





# ---- Crinier et al. 13-gene signature scoring ----

# Path to Excel file containing gene lists
FILE_SIGNATURES_PATH <- file.path(PATH_EXPERIMENT_REFERENCE, "Tables_Gene_Signatures", "Crinier_13_NKgenes.xlsx")

# ---- Function to read all sheets from Excel file ----
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if (!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

# ---- Load gene sets from all Excel sheets ----
MONITORED_Markers <- read_excel_allsheets(FILE_SIGNATURES_PATH)

# ---- Clean and truncate gene sets ----
# - Keep only genes found in the dataset
# - Limit to first 20 genes per sheet
MONITORED_Markers <- lapply(MONITORED_Markers, function(x) intersect(unlist(x), rownames(myObjectSeurat)))
MONITORED_Markers <- lapply(MONITORED_Markers, function(x) head(x, n = 20))

# ---- Score cells based on Crinier gene sets ----
myObjectSeurat <- AddModuleScore(
  myObjectSeurat,
  features = MONITORED_Markers,
  name = "Crinier_13genes_Score",
  pool = NULL,
  seed = 19
)

# ---- Plot the score distribution ----
p5 <- VlnPlot(
  myObjectSeurat,
  features = "Crinier_13genes_Score1",  # "1" added by Seurat
  y.max = MAX_CUTOFF_SCORE_MODULE,
  pt.size = 0
) & stat_summary(fun.data = data_summary, color = "black")

cat(' \n \n')
print(p5)
cat(' \n \n')



