## @knitr Azimuth_Annotation

# ---- Section header for RMarkdown output ----
cat(" \n \n")
cat("### BarPlots ")
cat(" \n \n")

# ---- Run Azimuth annotation ----
# Annotates query cells by projecting them onto a reference atlas 
myObjectSeurat <- RunAzimuth(myObjectSeurat, reference = Reference_Azimuth_To_Use)


# ---- Define color palettes ----

# Number of unique SingleR identities
nb_SingleR <- length(unique(myObjectSeurat$BlueprintEncodeData.fine))

# ---- SingleR palette ----
palette_SingleR <- hue_pal()(nb_SingleR)
names(palette_SingleR) <- unique(myObjectSeurat$BlueprintEncodeData.fine)
palette_SingleR["NK cells"] <- "#FF0000"  # Force NK cells to appear in red

# ---- Azimuth level 1 (main types) palette ----
palette_Azimuth1 <- hue_pal()(length(unique(myObjectSeurat$predicted.celltype.l1)))
names(palette_Azimuth1) <- unique(myObjectSeurat$predicted.celltype.l1)
palette_Azimuth1["NK"] <- "#FF0000"

# ---- Azimuth level 2 (fine types) palette ----
palette_Azimuth2 <- hue_pal()(length(unique(myObjectSeurat$predicted.celltype.l2)))
names(palette_Azimuth2) <- unique(myObjectSeurat$predicted.celltype.l2)
palette_Azimuth2["NK"] <- "#FF0000"


# ---- Generate barplots comparing cell type annotations across clusters ----

# Violin plot of NK gene module score
p1 <- VlnPlot(
  myObjectSeurat,
  features = "Crinier_13genes_Score1",
  y.max = 1.5,
  pt.size = 0
) & stat_summary(fun.data = data_summary, color = "black") & ggtitle("")

# Barplot: Cluster composition by SingleR labels
p2 <- ggplot(myObjectSeurat@meta.data, aes(x = seurat_clusters, fill = BlueprintEncodeData.fine)) +
  geom_bar(position = "fill") +
  theme(text = element_text(size = GGPLOT_ELEMENT_TEXT_SIZE)) +
  theme_bw() +
  scale_fill_manual(values = palette_SingleR)

# Barplot: Cluster composition by Azimuth L1 labels
p3 <- ggplot(myObjectSeurat@meta.data, aes(x = seurat_clusters, fill = predicted.celltype.l1)) +
  geom_bar(position = "fill") +
  theme(text = element_text(size = GGPLOT_ELEMENT_TEXT_SIZE)) +
  theme_bw() +
  scale_fill_manual(values = palette_Azimuth1)

# Barplot: Cluster composition by Azimuth L2 labels
p4 <- ggplot(myObjectSeurat@meta.data, aes(x = seurat_clusters, fill = predicted.celltype.l2)) +
  geom_bar(position = "fill") +
  theme(text = element_text(size = GGPLOT_ELEMENT_TEXT_SIZE)) +
  theme_bw() +
  scale_fill_manual(values = palette_Azimuth2)

# Combine the 4 barplots into one figure
figure_barplots <- plot_grid(
  p1, p2, p3, p4,
  labels = c("Crinier Score", "Single R", "Azimuth Main cell type", "Azimuth Precise Cell Type"),
  hjust = -0.75,
  ncol = 2, nrow = 2
)

cat(' \n \n')
print(figure_barplots)
cat(' \n \n')


# ---- UMAP in query space ----
cat(" \n \n")
cat("### UMAP ")
cat(" \n \n")

# UMAP colored by Seurat clusters
p1b <- DimPlot(myObjectSeurat, label = TRUE) & ggtitle("")

# UMAP colored by SingleR labels
p2b <- DimPlot(myObjectSeurat, group.by = "BlueprintEncodeData.fine", cols = palette_SingleR) & ggtitle("")

# UMAP colored by Azimuth L1 labels
p3b <- DimPlot(myObjectSeurat, group.by = "predicted.celltype.l1", cols = palette_Azimuth1) & ggtitle("")

# UMAP colored by Azimuth L2 labels
p4b <- DimPlot(myObjectSeurat, group.by = "predicted.celltype.l2", cols = palette_Azimuth2) & ggtitle("")

# Combine the 4 UMAPs into one figure
figure_UMAP <- plot_grid(
  p1b, p2b, p3b, p4b,
  labels = c("Clusters", "Single R", "Azimuth Main cell type", "Azimuth Precise Cell Type"),
  hjust = -0.75,
  ncol = 2, nrow = 2
)

cat(' \n \n')
print(figure_UMAP)
cat(' \n \n')


# ---- UMAP in Azimuth reference space ----
cat(" \n \n")
cat("### UMAP in reference UMAP ")
cat(" \n \n")

# UMAP in reference space (reduction = "ref.umap") colored by Seurat clusters
p1c <- DimPlot(myObjectSeurat, label = TRUE, reduction = "ref.umap") & ggtitle("")

# UMAP in reference space colored by SingleR labels
p2c <- DimPlot(myObjectSeurat, group.by = "BlueprintEncodeData.fine", reduction = "ref.umap", cols = palette_SingleR) & ggtitle("")

# UMAP in reference space colored by Azimuth L1 labels
p3c <- DimPlot(myObjectSeurat, group.by = "predicted.celltype.l1", reduction = "ref.umap", cols = palette_Azimuth1) & ggtitle("")

# UMAP in reference space colored by Azimuth L2 labels
p4c <- DimPlot(myObjectSeurat, group.by = "predicted.celltype.l2", reduction = "ref.umap", cols = palette_Azimuth2) & ggtitle("")

# Combine the 4 reference UMAPs into one figure
figure_UMAP_ref <- plot_grid(
  p1c, p2c, p3c, p4c,
  labels = c("Clusters", "Single R", "Azimuth Main cell type", "Azimuth Precise Cell Type"),
  hjust = -0.75,
  ncol = 2, nrow = 2
)

cat(' \n \n')
print(figure_UMAP_ref)
cat(' \n \n')





