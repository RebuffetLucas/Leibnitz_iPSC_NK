## @knitr correlation_meta_data_PC
# **Correlation of Principal Components with Metadata**
# This section computes the Spearman correlation between the principal components (PCs) and numeric metadata fields,
# followed by visualizing the correlation matrix.

cat(" \n \n")
cat("#### Correlation PCs with metadata ")
cat(" \n \n")

# Extract the principal component (PC) embeddings from the Seurat object
df <- Embeddings(myObjectSeurat)

# Combine the PC embeddings with the numeric columns of the metadata
# `select_if` is used to filter metadata columns that are numeric
df <- cbind(df, myObjectSeurat@meta.data[, colnames(select_if(myObjectSeurat@meta.data, is.numeric))])

# Compute the Spearman correlation matrix for the combined data
df <- cor(df, method = "spearman")

# Subset the correlation matrix to include only correlations between numeric metadata and the first 30 PCs
df <- df[
  colnames(select_if(myObjectSeurat@meta.data, is.numeric)),
  colnames(Embeddings(myObjectSeurat, reduction = "pca"))[1:30]
]

# Remove rows with only NA values from the correlation matrix
df <- df[rowSums(is.na(df)) != ncol(df), ]

# Visualize the correlation matrix using a correlation plot
# `tl.col` sets the text label color, `tl.srt` rotates the text labels, and `type = "full"` ensures a full correlation plot
corrplot(df, tl.col = "black", tl.srt = 45, type = "full", is.corr = TRUE)
