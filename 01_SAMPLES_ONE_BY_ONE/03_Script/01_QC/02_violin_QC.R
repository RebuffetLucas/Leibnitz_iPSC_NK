## @knitr violin_QC

# Print section header for RMarkdown/knitr output
cat(" \n \n")
cat("#### Violin QC ")
cat(" \n \n")

# Create a data frame with selected QC metrics from Seurat object
dfCont <- data.frame(
    nFeature_RNA = myObjectSeurat[["nFeature_RNA"]],      # Number of detected genes per cell
    nCount_RNA = myObjectSeurat[["nCount_RNA"]],          # Total UMI counts per cell
    percent.mito = myObjectSeurat[["percent.mito"]],      # Percent of reads mapping to mitochondrial genes
    orig.ident = myObjectSeurat[["orig.ident"]],          # Sample identifier
    percent.ribo = myObjectSeurat[["percent.ribo"]]       # Percent of reads mapping to ribosomal genes
)

# ---- Violin plot for number of features (genes detected) ----
p1 <- ggplot(dfCont, aes(x = orig.ident, y = nFeature_RNA, fill = orig.ident)) + 
    ggdist::stat_halfeye(           # Smoothed distribution plot
        adjust = .5, 
        width = .6, 
        .width = 0, 
        justification = -.3, 
        point_colour = NA
    ) + 
    geom_boxplot(                   # Add boxplot inside violin
        width = .15, 
        outlier.shape = NA
    ) +
    gghalves::geom_half_point(     # Add half-jittered points
        side = "l", 
        range_scale = .6, 
        alpha = .1,
        size = 0.1
    ) +
    scale_fill_manual(values = my_palette)  # Use custom color palette

# ---- Violin plot for total UMI counts ----
p2 <- ggplot(dfCont, aes(x = orig.ident, y = nCount_RNA, fill = orig.ident)) + 
    ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA) + 
    geom_boxplot(width = .15, outlier.shape = NA) +
    gghalves::geom_half_point(side = "l", range_scale = .6, alpha = .1, size = 0.1) +
    scale_fill_manual(values = my_palette)

# ---- Violin plot for percent mitochondrial content ----
p3 <- ggplot(dfCont, aes(x = orig.ident, y = percent.mito, fill = orig.ident)) + 
    ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA) + 
    geom_boxplot(width = .15, outlier.shape = NA) +
    gghalves::geom_half_point(side = "l", range_scale = .6, alpha = .1, size = 0.1) +
    scale_fill_manual(values = my_palette)

# ---- Violin plot for percent ribosomal content ----
p4 <- ggplot(dfCont, aes(x = orig.ident, y = percent.ribo, fill = orig.ident)) + 
    ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA) + 
    geom_boxplot(width = .15, outlier.shape = NA) +
    gghalves::geom_half_point(side = "l", range_scale = .6, alpha = .1, size = 0.1) +
    scale_fill_manual(values = my_palette)

# Combine and print all violin plots in one row
print(p1 + p2 + p3 + p4)


## @knitr FeatureScatter

cat(" \n \n")
cat("#### Feature Scatter ")
cat(" \n \n")

# ---- Scatter plot: nCount_RNA vs percent.mito ----
print(
    FeatureScatter(
        myObjectSeurat,
        feature1 = "nCount_RNA",
        feature2 = "percent.mito",
        group.by = "orig.ident",
        cols = my_palette,
        pt.size = PT_SIZE_FEATURE_SCATTER
    ) +
        geom_hline(yintercept = VISUAL_TRESHOLD_MITO, linetype = "dashed", color = "red")  # Horizontal QC threshold
)

# ---- Scatter plot: nCount_RNA vs nFeature_RNA ----
print(
    FeatureScatter(
        myObjectSeurat,
        feature1 = "nCount_RNA",
        feature2 = "nFeature_RNA",
        group.by = "orig.ident",
        cols = my_palette,
        pt.size = PT_SIZE_FEATURE_SCATTER
    )
)

# ---- Scatter plot: percent.mito vs percent.ribo ----
print(
    FeatureScatter(
        myObjectSeurat,
        feature1 = 'percent.mito',
        feature2 = 'percent.ribo',
        group.by = "orig.ident",
        cols = my_palette,
        pt.size = PT_SIZE_FEATURE_SCATTER
    ) +
        geom_vline(xintercept = VISUAL_TRESHOLD_MITO, linetype = "dashed", color = "red")  # Vertical mito QC threshold
)


############################################################
####### Same thing with SingleR naming as colour ###########
############################################################

## @knitr FeatureScatter2

# Define color palette for cell types from SingleR
Nb_pops = length(unique(myObjectSeurat@meta.data[["BlueprintEncodeData.fine"]]))    # Number of unique labels
mypalette2 = hue_pal()(Nb_pops)                                                     # Generate distinct colors
names(mypalette2) <- unique(myObjectSeurat@meta.data[["BlueprintEncodeData.fine"]]) # Name them
mypalette2["NK cells"] <- "#FF0000"                                                 # Assign red to NK cells

cat(" \n \n")
cat("#### Feature Scatter with SingleR Id")
cat(" \n \n")


# ---- Scatter plot: nCount_RNA vs nFeature_RNA, colored by cell type ----
print(
    FeatureScatter(
        myObjectSeurat,
        feature1 = "nCount_RNA",
        feature2 = "nFeature_RNA",
        group.by = "BlueprintEncodeData.fine",
        cols = mypalette2,
        pt.size = PT_SIZE_FEATURE_SCATTER
    ) +
        geom_hline(yintercept = VISUAL_TRESHOLD_nFEATURE, linetype = "dashed", color = "red") +
        geom_vline(xintercept = VISUAL_TRESHOLD_nCOUNT, linetype = "dashed", color = "red")
)

# ---- Scatter plot: percent.mito vs percent.ribo, colored by cell type ----
print(
    FeatureScatter(
        myObjectSeurat,
        feature1 = 'percent.mito',
        feature2 = 'percent.ribo',
        group.by = "BlueprintEncodeData.fine",
        cols = mypalette2,
        pt.size = PT_SIZE_FEATURE_SCATTER
    ) +
        geom_vline(xintercept = VISUAL_TRESHOLD_MITO, linetype = "dashed", color = "red")
)
