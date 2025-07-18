## @knitr Variable_Features

# Print section header for RMarkdown/knitr
cat(" \n \n")
cat("#### Variable Features ")
cat(" \n \n")

# ---- Normalize the data ----
# Applies log-normalization (counts per 10,000 + log1p transform)
myObjectSeurat <- NormalizeData(myObjectSeurat, normalization.method = "LogNormalize", scale.factor = SCALE_FACTOR)

# ---- Identify variable features ----
# Uses 'vst' method to find the most variable genes; nfeatures is user-defined (N_GENES_VARIABLES)
myObjectSeurat <- FindVariableFeatures(myObjectSeurat, selection.method = "vst", nfeatures = N_GENES_VARIABLES)

# ---- Extract the top 10 most variable genes ----
top10 <- head(VariableFeatures(myObjectSeurat), 10)

# ---- Plot variable features ----
plot1 <- VariableFeaturePlot(myObjectSeurat)            # Scatter plot of mean-variance relationship
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)  # Label top 10 most variable genes

print(plot2)  # Display plot with labeled genes





## @knitr scaling_and_pca

cat(" \n \n")
cat("#### Scaling and PCA ")
cat(" \n \n")

# ---- Extract gene names to scale all genes ----
all.genes <- rownames(myObjectSeurat)

# ---- Scale the data ----
# Centers and/or scales expression values (depending on DO_CENTER, DO_SCALE)
myObjectSeurat <- ScaleData(
  myObjectSeurat,
  features = all.genes,
  do.scale = DO_SCALE,     # Boolean: scale genes to unit variance
  do.center = DO_CENTER    # Boolean: center genes to mean zero
)

# ---- Run PCA ----
# Performs PCA using only the previously identified variable genes
myObjectSeurat <- RunPCA(myObjectSeurat, features = VariableFeatures(object = myObjectSeurat))

# ---- Visualize PCA result ----
print(DimPlot(myObjectSeurat, reduction = "pca"))  # Shows 2D PCA embedding colored by cluster or group



## @knitr UMAP

cat(" \n \n")
cat("#### UMAP ")
cat(" \n \n")

# ---- Construct k-nearest neighbor graph ----
# Uses PCA dimensions 1 to DIM_PCA (user-defined) as input
myObjectSeurat <- FindNeighbors(myObjectSeurat, dims = 1:DIM_PCA)

# ---- Cluster the cells ----
# Applies a graph-based clustering algorithm with user-defined resolution
myObjectSeurat <- FindClusters(myObjectSeurat, resolution = RESOLUTION, verbose = FALSE)

# ---- Run UMAP for visualization ----
# Uses PCA dimensions 1 to DIM_UMAP to compute UMAP embedding
myObjectSeurat <- RunUMAP(myObjectSeurat, dims = 1:DIM_UMAP)

# ---- Visualize UMAP ----
# Displays the UMAP embedding with clusters labeled
print(
  DimPlot(
    myObjectSeurat,
    reduction = "umap",
    pt.size = UMAP_PT_SIZE,     # Point size for cells
    label = TRUE,               # Show cluster labels
    label.size = LABEL_SIZE     # Font size of labels
  )
)


## @knitr markers

cat(" \n \n")
cat("#### Markers ")
cat(" \n \n")

# ---- Find marker genes for all clusters ----
# Identifies genes that are significantly more expressed in each cluster vs. all others
# only.pos = TRUE keeps only positive markers
# logfc.threshold = FIND_ALL_MARKERS_LOGFC is a user-defined log fold-change cutoff
markers <- FindAllMarkers(
  myObjectSeurat,
  only.pos = TRUE,
  logfc.threshold = FIND_ALL_MARKERS_LOGFC
)

# ---- Filter markers by adjusted p-value ----
# Keep only markers with FDR below threshold FDRCUTOFF
MarNKObject <- markers[markers$p_val_adj < FDRCUTOFF, ]

# ---- Prepare color palette ----
Nb_markers <- length(unique(myObjectSeurat@active.ident))  # Number of clusters
mypalette <- hue_pal()(Nb_markers)                         # Generate distinct colors

# ---- Display marker table with cluster-colored background ----
# Interactive table with download buttons and color-coded cluster column
print(
  htmltools::tagList(
    DT::datatable(
      MarNKObject,
      rownames = FALSE,
      extensions = 'Buttons',
      options = list(
        dom = 'Blfrtip',
        buttons = c('excel', "csv"),  # Enable export options
        fixedHeader = TRUE
      )
    ) %>%
      DT::formatStyle(
        'cluster',
        backgroundColor = DT::styleEqual(
          sort(unique(myObjectSeurat@active.ident)),
          mypalette[1:Nb_markers]
        )
      )
  )
)

