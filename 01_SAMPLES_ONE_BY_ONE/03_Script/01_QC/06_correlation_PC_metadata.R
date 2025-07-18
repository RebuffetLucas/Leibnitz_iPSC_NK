## @knitr correlation_meta_data_PC

# ---- Section header for RMarkdown/knitr ----
cat(" \n \n")
cat("#### Correlation PCs with metadata ")
cat(" \n \n")

# ---- Step 1: Retrieve PCA embeddings ----
# Extract principal component scores for each cell
df <- Embeddings(myObjectSeurat)  # This returns a matrix of PC scores (rows = cells, cols = PCs)

# ---- Step 2: Add numeric metadata to the dataframe ----
# Combine PC embeddings with numeric columns from metadata (e.g., nCount_RNA, percent.mito, etc.)
df <- cbind(
  df,
  myObjectSeurat@meta.data[, colnames(select_if(myObjectSeurat@meta.data, is.numeric))]
)

# ---- Step 3: Compute correlation matrix ----
# Compute Spearman correlation between metadata and PCs
df <- cor(df, method = "spearman")  # Rows and columns: metadata + PCs

# ---- Step 4: Subset correlation matrix ----
# Keep only correlations between:
# - rows: numeric metadata variables
# - columns: first 30 PCs
df <- df[
  colnames(select_if(myObjectSeurat@meta.data, is.numeric)),                      # Metadata rows
  colnames(Embeddings(myObjectSeurat, reduction = "pca"))[1:NUMBER_PCs_CORRELATION_PLOTS]                  # PC columns
]

# ---- Step 5: Remove rows with only NA values ----
df <- df[rowSums(is.na(df)) != ncol(df), ]

# ---- Step 6: Plot the correlation matrix ----
# Use corrplot for a visual representation (color-coded matrix)
corrplot(
  df,
  tl.col = "black",   # Color of text labels
  tl.srt = 45,        # Rotation angle of labels
  type = "full",      # Show full matrix
  is.corr = TRUE      # Matrix contains correlation values
)
