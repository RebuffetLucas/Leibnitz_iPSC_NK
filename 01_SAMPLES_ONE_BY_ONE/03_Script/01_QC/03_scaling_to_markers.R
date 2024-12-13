## @knitr Variable_Features

cat(" \n \n")
cat("#### Variable Features ")
cat(" \n \n")


myObjectSeurat <- NormalizeData(myObjectSeurat, normalization.method = "LogNormalize", scale.factor = 10000)
myObjectSeurat <- FindVariableFeatures(myObjectSeurat, selection.method = "vst", nfeatures = N_GENES_VARIABLES)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(myObjectSeurat), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(myObjectSeurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

print( plot2 )



## @knitr scaling_and_pca

cat(" \n \n")
cat("#### Scaling and PCA ")
cat(" \n \n")


all.genes <- rownames(myObjectSeurat)
myObjectSeurat <- ScaleData(myObjectSeurat, features = all.genes,do.scale = DO_SCALE,do.center = DO_CENTER)
myObjectSeurat <- RunPCA(myObjectSeurat, features = VariableFeatures(object = myObjectSeurat))
print(DimPlot(myObjectSeurat, reduction = "pca"))


## @knitr UMAP

cat(" \n \n")
cat("#### UMAP ")
cat(" \n \n")


myObjectSeurat <- FindNeighbors(myObjectSeurat, dims = 1:DIM_PCA)
myObjectSeurat <- FindClusters(myObjectSeurat, resolution = RESOLUTION, verbose = FALSE)
myObjectSeurat <- RunUMAP(myObjectSeurat, dims = 1:DIM_UMAP)

print(DimPlot(myObjectSeurat, reduction = "umap", pt.size = UMAP_PT_SIZE, label = TRUE, label.size = LABEL_SIZE))


## @knitr markers


cat(" \n \n")
cat("#### Markers ")
cat(" \n \n")

markers <- FindAllMarkers(myObjectSeurat , only.pos = TRUE,logfc.threshold = FIND_ALL_MARKERS_LOGFC)
MarNKObject=(markers[markers$p_val_adj<FDRCUTOFF,])
	
Nb_markers=length(unique(myObjectSeurat@active.ident))
mypalette= hue_pal()(Nb_markers)

 
print( htmltools::tagList(DT::datatable(MarNKObject,rownames = FALSE,extensions = 'Buttons', 
	    options = list(dom = 'Blfrtip', 
		           buttons = c('excel', "csv"), fixedHeader = TRUE)
	    )%>% 
	  DT::formatStyle(
	  'cluster',
	  backgroundColor = DT::styleEqual(sort(unique(myObjectSeurat@active.ident)),  mypalette[1:Nb_markers])
	)))
