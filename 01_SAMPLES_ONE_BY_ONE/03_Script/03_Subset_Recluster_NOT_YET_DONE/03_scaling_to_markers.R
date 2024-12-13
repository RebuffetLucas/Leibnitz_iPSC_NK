## @knitr Variable_Features
# **Identification of Highly Variable Features**
# This section identifies and visualizes the most variable features in the dataset.

cat(" \n \n")
cat("#### Variable Features ")
cat(" \n \n")

# Normalize the data using log-normalization
myObjectSeurat <- NormalizeData(myObjectSeurat, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify variable features using the "vst" method and retain the top 'n' features
myObjectSeurat <- FindVariableFeatures(myObjectSeurat, selection.method = "vst", nfeatures = N_GENES_VARIABLES)

# Extract the top 10 most variable features
top10 <- head(VariableFeatures(myObjectSeurat), 10)

# Plot variable features and highlight the top 10 with labels
plot1 <- VariableFeaturePlot(myObjectSeurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

# Display the annotated plot
print(plot2)

## @knitr scaling_and_pca
# **Scaling and Principal Component Analysis (PCA)**
# This section scales the data and performs PCA to reduce dimensionality.

cat(" \n \n")
cat("#### Scaling and PCA ")
cat(" \n \n")

# Retrieve all gene names
all.genes <- rownames(myObjectSeurat)

# Scale the data by centering and scaling each feature
myObjectSeurat <- ScaleData(myObjectSeurat, features = all.genes, do.scale = DO_SCALE, do.center = DO_CENTER)

# Run PCA on the scaled data using the identified variable features
myObjectSeurat <- RunPCA(myObjectSeurat, features = VariableFeatures(object = myObjectSeurat))

# Display the PCA plot
print(DimPlot(myObjectSeurat, reduction = "pca"))

## @knitr UMAP
# **UMAP Dimensionality Reduction and Clustering**
# This section performs UMAP for visualization and identifies clusters.

cat(" \n \n")
cat("#### UMAP ")
cat(" \n \n")

# Find neighbors based on the top PCA dimensions
myObjectSeurat <- FindNeighbors(myObjectSeurat, dims = 1:DIM_PCA)

# Identify clusters using the Louvain algorithm with the specified resolution
myObjectSeurat <- FindClusters(myObjectSeurat, resolution = RESOLUTION, verbose = FALSE)

# Run UMAP for visualization
myObjectSeurat <- RunUMAP(myObjectSeurat, dims = 1:DIM_UMAP)

# Display the UMAP plot with cluster labels
print(DimPlot(myObjectSeurat, reduction = "umap", pt.size = UMAP_PT_SIZE, label = TRUE, label.size = LABEL_SIZE))

## @knitr markers
# **Identification and Visualization of Cluster-Specific Markers**
# This section identifies markers for each cluster and creates an interactive table.

cat(" \n \n")
cat("#### Markers ")
cat(" \n \n")

# Identify markers for all clusters, filtering for positive markers with a log fold change threshold
markers <- FindAllMarkers(myObjectSeurat, only.pos = TRUE, logfc.threshold = FIND_ALL_MARKERS_LOGFC)

# Filter markers based on adjusted p-value (FDR cutoff)
MarNKObject <- markers[markers$p_val_adj < FDRCUTOFF,]

# Retrieve the number of unique clusters
Nb_markers <- length(unique(myObjectSeurat@active.ident))

# Define a color palette for clusters
mypalette <- hue_pal()(Nb_markers)

# Create an interactive table displaying the filtered markers with cluster-specific colors
print(htmltools::tagList(DT::datatable(
  MarNKObject, 
  rownames = FALSE, 
  extensions = 'Buttons',
  options = list(dom = 'Blfrtip', buttons = c('excel', "csv"), fixedHeader = TRUE)
) %>% 
  DT::formatStyle(
    'cluster',
    backgroundColor = DT::styleEqual(sort(unique(myObjectSeurat@active.ident)), mypalette[1:Nb_markers])
  )))
