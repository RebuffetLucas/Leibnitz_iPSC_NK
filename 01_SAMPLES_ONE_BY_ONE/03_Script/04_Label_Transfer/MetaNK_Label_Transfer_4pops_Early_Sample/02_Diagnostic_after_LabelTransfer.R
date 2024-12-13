## @knitr DiagnosticAfterLabelTransfer
# **Diagnostic After Label Transfer**
# This script performs a series of diagnostic visualizations and quality assessments 
# to evaluate the accuracy and quality of label transfer and the resulting predictions.

# Display tabbed structure for diagnostics
cat("# Diagnostic After Label Transfer {.tabset .tabset-fade} \n\n")

# **Global UMAP Plots**
# Visualize UMAP representations of the reference dataset and query dataset.
cat("## Diag heterogeneity quality Data {.tabset .tabset-fade} \n\n")
cat("### Data and Meta NK data raw \n\n")

p01 = DimPlot(PBMC_Meta_NK_V2, cols = palette2)  # UMAP of reference dataset with cluster colors
p02 = DimPlot(PBMC_Blood_All, group.by = QUERY_CLUSTER_COLNAME)  # Query dataset by cluster
p03 = DimPlot(PBMC_Blood_All, group.by = QUERY_ORIG_IDENT_COLNAME)  # Query dataset by original identity

cat("\n\n")
print(p01)
cat("\n\n")
cat("\n\n")
print(p02 + p03)
cat("\n\n")

# **Diagnostic of Counts and Features**
# Explore RNA counts and detected features across clusters or original identities.
cat("### Diag data Counts and Features \n\n")

PBMC_Blood_All$nCount = colSums(x = PBMC_Blood_All, slot = "counts")  # Total RNA counts
PBMC_Blood_All$nFeature = colSums(x = GetAssayData(object = PBMC_Blood_All, slot = "counts") > 0)  # Number of features

p04 = VlnPlot(PBMC_Blood_All, features = "nCount", group.by = QUERY_CLUSTER_COLNAME) & stat_summary(fun.data = data_summary, color = "black")
p05 = VlnPlot(PBMC_Blood_All, features = "nCount", group.by = QUERY_ORIG_IDENT_COLNAME) & stat_summary(fun.data = data_summary, color = "black")

cat("\n\n")
print(p04)
cat("\n\n")
cat("\n\n")
print(p05)
cat("\n\n")

# **Visualization of Reference**
# UMAP visualizations of the integrated reference dataset.
cat("## Vizualization of the reference {.tabset .tabset-fade} \n\n")


p1 <- DimPlot(data_integrated, reduction = "umap", group.by = "Dataset")
p2 <- DimPlot(data_integrated, reduction = "umap", group.by = "FirstClust", label = TRUE, 
              repel = TRUE, cols = palette2) + NoLegend()

cat("\n\n")
print( p1 + p2 )
cat("\n\n")

p4 <- DimPlot(data_integrated, reduction = "umap", split.by = "Dataset", group.by = "FirstClust", label = TRUE, 
              repel = TRUE, cols = palette2) + NoLegend()

cat("\n\n")
print(p4)
cat("\n\n")


# **Assessing Quality of Label Transfer**
# Analyze prediction reliability and marker consistency across predicted populations.
cat("## Assessing quality of Label Transfer {.tabset .tabset-fade} \n\n")

cat("### Look at the common markers for the different populations {.tabset .tabset-fade} \n\n")

data.query = SetIdent(data.query, value = "predicted.id")

cat("#### Main pops based on V2 \n\n")


#Main pops
Markers_NK_123 = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/TOP20_MainClusters.rds")
  
# **Violin Plot of Prediction Scores**
# Remove small clusters to improve visualization.
if (DO_REMOVE_VERY_LITTLE_IDENT_SUBSET==TRUE){
  # Get the counts of each predicted.id
  predicted_counts <- table(data.query$predicted.id)
  # Retrieve the names with counts less than TRESHOLD_NEW_ID_REMOVING
  predicted_less_than_treshold <- names(predicted_counts[predicted_counts < TRESHOLD_NEW_ID_REMOVING])
  data.query2 = subset(data.query, subset= predicted.id %in% predicted_less_than_treshold , invert=TRUE)
}else{data.query2 = data.query
}

# Add module scores for markers and plot violin plots of prediction scores.
for (i in names(Markers_NK_123)){
  data.query2 = AddModuleScore(data.query2, features = as.list(Markers_NK_123[[i]]), pool= NULL ,name= i , seed=19)
}

p_Vln1 = VlnPlot(data.query2, features = paste0(names(Markers_NK_123),"1") , group.by = "predicted.id", cols= palette2, pt.size = 0)  & stat_summary(fun.data=data_summary,color="black")

cat("\n\n")
print(p_Vln1)
cat("\n\n")


# **Dot Plots for Markers**
# Visualize marker expression for each predicted population.
p5 = DotPlot(data.query2, features = Markers_NK_123$NK1, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK1 genes") + theme(axis.text = element_text(size = 20)) 
p6 = DotPlot(data.query2, features = Markers_NK_123$NK2, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK2 genes") + theme(axis.text = element_text(size = 20)) 
p7 = DotPlot(data.query2, features = Markers_NK_123$NK3, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK3 genes") + theme(axis.text = element_text(size = 20)) 

cat("\n\n")
print(p5)
cat("\n\n")
print(p6)
cat("\n\n")
print(p7)
cat("\n\n")

# **Reliability of Predictions**
# Visualize prediction scores for each population and across different clusters.

cat("#### Main pops based on CITEseq \n\n")

#Main pops
Markers_NK_123_CITEseq = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/DEG/All_Markers_CITEseq3clusters.rds")

Markers_NK_123_CITEseq %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  filter(!grepl("RPL|RPS|MT-", gene)) %>% #We use RPS and RPL for scoring but nor for plotind DotPlots
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_TOP_SCORING, wt = avg_log2FC) -> top_All

top_All %>% filter(cluster== "NK_1") -> Top_NK1
top_All %>% filter(cluster== "NK_2") -> Top_NK2
top_All %>% filter(cluster== "NK_3") -> Top_NK3


#ViolinPlot
#Remove the very little identified clusters that could mislead dataviz in the scaled plots
if (DO_REMOVE_VERY_LITTLE_IDENT_SUBSET==TRUE){
  # Get the counts of each predicted.id
  predicted_counts <- table(data.query$predicted.id)
  # Retrieve the names with counts less than TRESHOLD_NEW_ID_REMOVING
  predicted_less_than_treshold <- names(predicted_counts[predicted_counts < TRESHOLD_NEW_ID_REMOVING])
  data.query2 = subset(data.query, subset= predicted.id %in% predicted_less_than_treshold , invert=TRUE)
}else{data.query2 = data.query
}

data.query2 = AddModuleScore(data.query2, features = list(Top_NK1$gene), pool= NULL ,name= "NK1" , seed=19)
data.query2 = AddModuleScore(data.query2, features = list(Top_NK2$gene), pool= NULL ,name= "NK2" , seed=19)
data.query2 = AddModuleScore(data.query2, features = list(Top_NK3$gene), pool= NULL ,name= "NK3" , seed=19)


p_Vln2 = VlnPlot(data.query2, features = c("NK11", "NK21", "NK31") , group.by = "predicted.id", cols= palette2, pt.size = 0)  & stat_summary(fun.data=data_summary,color="black")

cat("\n\n")
print(p_Vln2)
cat("\n\n")

p8 = DotPlot(data.query2, features = Top_NK1$gene, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK1 genes") + theme(axis.text = element_text(size = 20)) 
p9 = DotPlot(data.query2, features = Top_NK2$gene, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK2 genes") + theme(axis.text = element_text(size = 20)) 
p10 = DotPlot(data.query2, features = Top_NK3$gene, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK3 genes") + theme(axis.text = element_text(size = 20)) 

cat("\n\n")
print(p8)
cat("\n\n")
print(p9)
cat("\n\n")
print(p10)
cat("\n\n")


cat("#### Secondary pops \n\n")
  #Secondary pops
Markers_NK_6pops = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/TOP11_6_Clusters.rds")

#Remove the very little identified clusters that could mislead dataviz in the scaled plots
if (DO_REMOVE_VERY_LITTLE_IDENT_SUBSET==TRUE){
  # Get the counts of each predicted.id
  predicted_counts <- table(data.query$predicted.id)
  # Retrieve the names with counts less than TRESHOLD_NEW_ID_REMOVING
  predicted_less_than_treshold <- names(predicted_counts[predicted_counts < TRESHOLD_NEW_ID_REMOVING])
  data.query2 = subset(data.query, subset= predicted.id %in% predicted_less_than_treshold , invert=TRUE)
}else{data.query2 = data.query
}

for (i in names(Markers_NK_6pops)){
  data.query2 = AddModuleScore(data.query2, features = as.list(Markers_NK_6pops[[i]]), pool= NULL ,name= i , seed=19)
}

p_Vln3 = VlnPlot(data.query2, features = paste0(names(Markers_NK_6pops),"1") , group.by = "predicted.id", cols= palette2, pt.size = 0)  & stat_summary(fun.data=data_summary,color="black") & stat_summary(fun.data=data_summary,color="black")

cat("\n\n")
print(p_Vln3)
cat("\n\n")

p11 = DotPlot(data.query2, features = Markers_NK_6pops$NK1A , cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK1A genes") + theme(axis.text = element_text(size = 20)) 
p12 = DotPlot(data.query2, features = Markers_NK_6pops$NK1B, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK1B genes") + theme(axis.text = element_text(size = 20)) 
p13 = DotPlot(data.query2, features = Markers_NK_6pops$NK1C, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK1C genes") + theme(axis.text = element_text(size = 20)) 
p14 = DotPlot(data.query2, features = Markers_NK_6pops$NKint, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NKint genes") + theme(axis.text = element_text(size = 20)) 
p15 = DotPlot(data.query2, features = Markers_NK_6pops$NK2, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK2 genes") + theme(axis.text = element_text(size = 20)) 
p16 = DotPlot(data.query2, features = Markers_NK_6pops$NK3, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK3 genes") + theme(axis.text = element_text(size = 20)) 

cat("\n\n")
print(p11)
cat("\n\n")
print(p12)
cat("\n\n")
print(p13)
cat("\n\n")
print(p14)
cat("\n\n")
print(p15)
cat("\n\n")
print(p16)
cat("\n\n")


# **Reliability of Predictions**
# Visualize prediction scores for each population and across different clusters.
cat("### Reliability of the prediction across predicted id \n\n")

p17 = VlnPlot(data.query2 , features= "prediction.score.max", group.by = "predicted.id" , cols = palette2, pt.size = 0) & stat_summary(fun.data=data_summary,color="black")

data.query2$predicted.id = as.factor(data.query2$predicted.id)




for(pop_predicted in levels(data.query2$predicted.id)){
  data.query2_subset = subset(data.query2, subset = predicted.id == pop_predicted)
  Vlnp_loop= VlnPlot(data.query2_subset , features= "prediction.score.max", group.by = QUERY_CLUSTER_COLNAME , pt.size = 0) & stat_summary(fun.data=data_summary,color="black") & ggtitle(pop_predicted)
  Ridge_loop= RidgePlot(data.query2_subset , features= "prediction.score.max", group.by = QUERY_CLUSTER_COLNAME )  & ggtitle(pop_predicted)
  
  cat("\n\n")
  print(Vlnp_loop)
  cat("\n\n")
  
  cat("\n\n")
  print(Ridge_loop)
  cat("\n\n")
  
  Vlnp_loop2= VlnPlot(data.query2_subset , features= "prediction.score.max", group.by = QUERY_ORIG_IDENT_COLNAME , pt.size = 0) & stat_summary(fun.data=data_summary,color="black") & ggtitle(pop_predicted)
  Ridge_loop2= RidgePlot(data.query2_subset , features= "prediction.score.max", group.by = QUERY_ORIG_IDENT_COLNAME )  & ggtitle(pop_predicted)
  
  cat("\n\n")
  print(Vlnp_loop2)
  cat("\n\n")
  
  cat("\n\n")
  print(Ridge_loop2)
  cat("\n\n")

}


cat("\n\n")
print(p17)
cat("\n\n")

p17b = VlnPlot(data.query2 , features= "prediction.score.max", group.by = QUERY_CLUSTER_COLNAME, pt.size = 0)  

cat("\n\n")
print(p17b)
cat("\n\n")

p17c = VlnPlot(data.query2 , features= "prediction.score.max", group.by = "predicted.id", split.by = QUERY_CLUSTER_COLNAME)  

cat("\n\n")
print(p17c)
cat("\n\n")

p17d = VlnPlot(data.query2 , features= "prediction.score.max", group.by = QUERY_CLUSTER_COLNAME, split.by = "predicted.id" , cols = c(palette2[1],palette2[3] ,palette2[4],palette2[2]))  
p17d

cat("\n\n")
print(p17d)
cat("\n\n")





cat("### Look at the spontaneous markers for the different populations predicted \n\n")
cat("#####" ,"DotPlot spontaneous markers","\n")
if (DO_REMOVE_VERY_LITTLE_IDENT_SUBSET==TRUE){
  # Get the counts of each predicted.id
  predicted_counts <- table(data.query$predicted.id)
  # Retrieve the names with counts less than TRESHOLD_NEW_ID_REMOVING
  predicted_less_than_treshold <- names(predicted_counts[predicted_counts < TRESHOLD_NEW_ID_REMOVING])
  data.query2 = subset(data.query, subset= predicted.id %in% predicted_less_than_treshold , invert=TRUE)
}else{data.query2 = data.query
}

data.query2 = SetIdent(data.query2, value = "predicted.id")

All_Markers = FindAllMarkers(data.query2 , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

#Look at markers
All_Markers %>%
  group_by(cluster) %>%
  filter(!grepl("RPL|RPS|MT-", gene)) %>%
  top_n(n = FINDMARKERS_SHOWTOP, wt = avg_log2FC) -> top10


p18 = DotPlot(data.query2, features =  unique(top10$gene) , cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + theme(axis.text = element_text(size = 20)) 

cat("\n\n")
print(p18)
cat("\n\n")

#interactive markers table
#Make the interactive table
DEG_sample <- All_Markers
DEG_sample_filtered =(DEG_sample[DEG_sample$p_val_adj<FDRCUTOFF,])
DEG_sample_filtered <- DEG_sample_filtered[order(DEG_sample_filtered$avg_log2FC ,decreasing = T),]

cat("#####" ,"Markers Table","\n")

Nb_markers=length(unique(data.query2@active.ident))
mypalette= hue_pal()(Nb_markers)
mypalette=palette2 #Vérifier si ça marche ou si ça fait bugger
print( htmltools::tagList(DT::datatable(DEG_sample_filtered,rownames = FALSE,extensions = 'Buttons', 
                                        options = list(dom = 'Blfrtip', 
                                                       buttons = c('excel', "csv"), fixedHeader = TRUE)
)%>% 
  DT::formatStyle(
    'cluster',
    backgroundColor = DT::styleEqual(sort(unique(data.query2@active.ident)),  mypalette[1:Nb_markers])
  )))


for(clust in sort(unique(Idents(data.query2)))){
  cat("\n\n")
  cat("##### Cluster ", clust,  "{.tabset .tabset-fade}\n\n")
  cat("\n\n")
  DEG_clust <- DEG_sample_filtered[DEG_sample_filtered$cluster %in% clust,]
  
  for(gene in head(DEG_clust$gene)){
    cat("\n\n")
    cat("######", gene)
    cat("\n\n")
    print(VlnPlot(data.query2, group.by = "predicted.id", features = gene, pt.size = 0))
    cat("\n\n")
  }
}



cat("## Look at the results of label transfer {.tabset .tabset-fade} \n\n")

cat("\n\n")
p22 = ggplot(data.query2@meta.data, aes(x=SecondClust, fill= predicted.id)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= palette2) + ggtitle("Data Only") + theme(text = element_text(size = 20))

cat("\n\n")
#print(p22)

cat("\n\n")
p21 = ggplot(data.query2@meta.data, aes_string(x=QUERY_CLUSTER_COLNAME, fill= "predicted.id")) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= palette2) + ggtitle("Data Only") + theme(text = element_text(size = 20))

cat("\n\n")
print(p21)


cat("\n\n")
p23 = ggplot(PBMC_Meta_NK_V2@meta.data, aes(x=orig.ident, fill= FirstClust)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90)) + scale_fill_manual(values= palette2) + ggtitle("Meta NK per Sample") + theme(text = element_text(size = 20))
print(p23)
cat("\n\n")

cat("\n\n")
p24 = ggplot(PBMC_Meta_NK_V2@meta.data, aes(x=Dataset, fill= FirstClust)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90)) + scale_fill_manual(values= palette2)  + ggtitle("Meta NK per Dataset") + theme(text = element_text(size = 20))
print(p24)
cat("\n\n")

