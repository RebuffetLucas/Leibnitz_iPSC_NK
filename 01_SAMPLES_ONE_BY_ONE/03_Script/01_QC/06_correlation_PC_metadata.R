## @knitr correlation_meta_data_PC

cat(" \n \n")
cat("#### Correlation PCs with metadata ")
cat(" \n \n")


df<-Embeddings(myObjectSeurat)
df <- cbind(df, myObjectSeurat@meta.data[,colnames(select_if(myObjectSeurat@meta.data, is.numeric))])
  
df <- cor(df, method = "spearman")
df <- df[colnames(select_if(myObjectSeurat@meta.data, is.numeric)),colnames(Embeddings(myObjectSeurat, reduction = "pca"))[1:30]]
df <- df[rowSums(is.na(df)) != ncol(df), ]
  
corrplot(df, tl.col = "black", tl.srt = 45, type = "full", is.corr = TRUE)

  