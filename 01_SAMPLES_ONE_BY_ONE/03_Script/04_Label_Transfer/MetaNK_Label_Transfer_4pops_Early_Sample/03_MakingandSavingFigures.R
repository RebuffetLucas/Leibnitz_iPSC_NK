## @knitr MakingandSavingFigures

cat("## Look at the results of label transfer in the initial UMAP {.tabset .tabset-fade} \n\n")
cat("### UMAP and Barplot \n\n")
#Reinject the prdicted.id in the initial object
predicted_id_to_add = readRDS(file.path(PATH_ANALYSIS_OUTPUT, "predicted_id_4pops.rds"))
names(predicted_id_to_add) <- sub("^Query_", "", names(predicted_id_to_add))

# Ensure the row names match
common_names <- intersect(rownames(PBMC_Blood_All@meta.data), names(predicted_id_to_add))

# Add the "predicted.id" column
PBMC_Blood_All@meta.data$predicted.id <- NA  # Initialize the column with NA
PBMC_Blood_All@meta.data[common_names, "predicted.id"] <- predicted_id_to_add[common_names]
p31 = DimPlot(PBMC_Blood_All, group.by = "predicted.id", pt.size = 0.6)

cat("\n\n")
print(p31)
cat("\n\n")


#p27 = ggplot(PBMC_Blood_All@meta.data, aes_string(x="RNA_snn_res.0.2", fill= "predicted.id")) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= palette2) + ggtitle("Data Only") + theme_minimal() + theme(text = element_text(size = 20))


p21 = ggplot(PBMC_Blood_All@meta.data, aes_string(x=QUERY_CLUSTER_COLNAME, fill= "predicted.id")) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= palette2) + ggtitle("Data Only") + theme_minimal() + theme(text = element_text(size = 20))

p21 <- ggplot(PBMC_Blood_All@meta.data, 
              aes_string(x = QUERY_CLUSTER_COLNAME, fill = "predicted.id")) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = palette2) +
  ggtitle("Data Only") +
  theme_minimal() +  # Use a minimal theme for a cleaner look
  theme(
    text = element_text(size = 20),       # General text size
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  # Rotate x-axis labels
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center-align and bold title
    legend.title = element_text(face = "bold"),  # Bold legend title
    legend.position = "right"  # Place legend on the right
  )


cat("\n\n")
print(p21)
cat("\n\n")

#Save the figures


pdf(file.path(PATH_ANALYSIS_OUTPUT_FIGURES, "UMAP_Predicted_Id.pdf"),  width = 10, height = 6)
p31
dev.off()

pdf(file.path(PATH_ANALYSIS_OUTPUT_FIGURES, "Barplot_Predicted_Id.pdf"),  width = 10, height = 6)
p21
dev.off()

DimPlot(PBMC_Blood_All , group.by = "RNA_snn_res.0.2")

#pdf(file.path(PATH_ANALYSIS_OUTPUT_FIGURES, "Barplot_oldclust_Predicted_Id.pdf"),  width = 10, height = 6)
#p27
#dev.off()


cat("### Alluvial \n\n")

#Plot and save the alluvial plot
# Load necessary libraries
library(ggalluvial)

if (DO_REMOVE_VERY_LITTLE_IDENT_SUBSET==TRUE){
  # Get the counts of each predicted.id
  predicted_counts <- table(PBMC_Blood_All$predicted.id)
  # Retrieve the names with counts less than TRESHOLD_NEW_ID_REMOVING
  predicted_less_than_treshold <- names(predicted_counts[predicted_counts < TRESHOLD_NEW_ID_REMOVING])
  PBMC_Blood_All2 = subset(PBMC_Blood_All, subset= predicted.id %in% predicted_less_than_treshold , invert=TRUE)
}else{PBMC_Blood_All2 = PBMC_Blood_All
}



# Assuming your data frame is called `meta_data`
PBMC_Blood_All2$predicted.id = as.factor(PBMC_Blood_All2$predicted.id )
meta_data <- PBMC_Blood_All2@meta.data

# Prepare the data for the alluvial plot
alluvial_data <- meta_data %>%
  dplyr::select(renum, predicted.id) %>%
  dplyr::count(renum, predicted.id)

# Create the alluvial plot
p32 = ggplot(alluvial_data, aes(
  axis1 = renum,
  axis2 = predicted.id,
  y = n
)) +
  geom_alluvium(aes(fill = renum)) +
  geom_stratum(width = 1/12, fill = "white") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("celltype", "predicted.id"), expand = c(0, 0)) +
  labs(
    title = "Alluvial Plot of celltype vs predicted.id",
    x = "Categories",
    y = "Frequency"
  ) +
  theme_minimal()

cat("\n\n")
print(p32)
cat("\n\n")


pdf(file.path(PATH_ANALYSIS_OUTPUT_FIGURES, "RidgePlot_Predicted_Id.pdf"),  width = 6, height = 5)
p32
dev.off()


#Plot and save the alluvial plot
cat("\n\n")
print(p32)
cat("\n\n")


pdf(file.path(PATH_ANALYSIS_OUTPUT_FIGURES, "RidgePlot_Predicted_Id.pdf"),  width = 6, height = 5)
p32
dev.off()


cat("### Scoring on UMAP \n\n")

Markers_NK_123 = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/TOP20_MainClusters.rds")

for (i in names(Markers_NK_123)){
  PBMC_Blood_All = AddModuleScore(PBMC_Blood_All, features = as.list(Markers_NK_123[[i]]), pool= NULL ,name= i , seed=19)
}

#Markers_NK_6pops = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/TOP11_6_Clusters.rds")
#names(Markers_NK_6pops[["NKint"]]) = "NKint"
#PBMC_Blood_All = AddModuleScore(PBMC_Blood_All, features = Markers_NK_6pops[["NKint"]], pool= NULL, name= "NKint", seed= 19)

#Ploting
features_to_plot <- paste0(names(Markers_NK_123), "1")
#features_to_plot = append(features_to_plot, "NKint1")
plot_titles <- c("Score NK1", "Score NK3" , "Score NK2")

#plot_titles <- c("Score NK1", "Score NK3" , "Score NK2", "Score NKint")

# Generate the FeaturePlots with titles
plots <- lapply(seq_along(features_to_plot), function(i) {
  FeaturePlot(PBMC_Blood_All, features = features_to_plot[i], pt.size = 0.6) &
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) &
    ggtitle(plot_titles[i])
})


# Combine the plots
combined_plot <- wrap_plots(plots)

# Display the combined plot
print(combined_plot)

pdf(file.path(PATH_ANALYSIS_OUTPUT_FIGURES, "NK123,scores_UMAP.pdf"),  width = 20, height = 4)
print(combined_plot)
dev.off()



#Remove the very little identified clusters that could mislead dataviz in the scaled plots
if (DO_REMOVE_VERY_LITTLE_IDENT_SUBSET==TRUE){
  # Get the counts of each predicted.id
  predicted_counts <- table(data.query$predicted.id)
  # Retrieve the names with counts less than TRESHOLD_NEW_ID_REMOVING
  predicted_less_than_treshold <- names(predicted_counts[predicted_counts < TRESHOLD_NEW_ID_REMOVING])
  data.query2 = subset(data.query, subset= predicted.id %in% predicted_less_than_treshold , invert=TRUE)
}else{data.query2 = data.query
}


cat("### Heatmap Overlap and Jaccard \n\n")


#Extract the genes
All_Markers_ref = readRDS( "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/signatures/FindAllMarkers_4pops_Output.rds")

#Look at markers
All_Markers_ref %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  #filter(!grepl("RPL|RPS|MT-", gene)) %>% #We use RPS and RPL for scoring but nor for plotind DotPlots
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_GENES_JACKARD, wt = avg_log2FC) -> top_ref

table(top_ref$cluster)

#Extract the gene lists
#PBMC_Blood_All_save = PBMC_Blood_All
#PBMC_Blood_All = subset(PBMC_Blood_All, features= genes_to_keep) #Keep only the genes that are shared
#PBMC_Blood_All = SetIdent(PBMC_Blood_All, value= "renum")

All_Markers_query = FindAllMarkers(PBMC_Blood_All , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

#Look at markers
All_Markers_query %>%
  group_by(cluster) %>%
  #filter(!grepl("RPL|RPS|MT-", gene)) %>%
  top_n(n = FINDMARKERS_SHOWTOP, wt = avg_log2FC) -> top_query

All_Markers_query %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  #filter(!grepl("RPL|RPS|MT-", gene)) %>% #We use RPS and RPL for scoring but nor for plotind DotPlots
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_GENES_JACKARD, wt = avg_log2FC) -> top_query

table(top_query$cluster)



#Jackard and overlap Index for 4 clusters
# Split the data frame into a list of data frames by cluster
list_per_cluster_ref <- split(top_ref, top_ref$cluster)
qlist_per_cluster_query <- split(top_query, top_query$cluster)
# Extract the gene column from each data frame in the list
genes_per_cluster_ref <- lapply(list_per_cluster_ref, function(df) df$gene)
genes_per_cluster_query <- lapply(qlist_per_cluster_query, function(df) df$gene)


# Function to calculate Jaccard Index
jaccard_index <- function(set1, set2) {
  if (length(set1) == 0 || length(set2) == 0) {
    return(0)
  }
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(intersection / union)
}

# Assuming genes_per_cluster_ref and genes_per_cluster_query are your lists of lists
# Calculate Jaccard Index for each combination of clusters between the two lists of lists
jaccard_results <- matrix(0, nrow = length(genes_per_cluster_ref), ncol = length(genes_per_cluster_query))
rownames(jaccard_results) <- names(genes_per_cluster_ref)
colnames(jaccard_results) <- names(genes_per_cluster_query)

for (i in seq_along(genes_per_cluster_ref)) {
  for (j in seq_along(genes_per_cluster_query)) {
    jaccard_results[i, j] <- jaccard_index(genes_per_cluster_ref[[i]], genes_per_cluster_query[[j]])
  }
}

# Convert the matrix to a long format suitable for ggplot
jaccard_data_long <- as.data.frame(as.table(jaccard_results))

# Renaming the columns for clarity
colnames(jaccard_data_long) <- c("Cluster_ref", "Cluster_query", "OverlapIndex")

# Create the heatmap
heat_jackard_ref_query =    ggplot(jaccard_data_long, aes(x = Cluster_ref, y = Cluster_query, fill = OverlapIndex)) +
  geom_tile(color="white") + # Use geom_tile() for heatmap squares
  scale_fill_viridis() + 
  theme_minimal() + # Use a minimal theme for cleaner appearance
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels for better readability
  labs(fill = "jaccard\nIndex", x = NULL, y = NULL, title = "jaccard Index Heatmap") # Set labels and title

print(heat_jackard_ref_query)

### OVERLAP ###
# Function to calculate Overlap Index
overlap_index <- function(set1, set2) {
  if (length(set1) == 0 || length(set2) == 0) {
    return(0)
  }
  intersection_size <- length(intersect(set1, set2))
  min_size <- min(length(set1), length(set2))
  return(intersection_size / min_size)
}

# Assuming genes_per_cluster_ref and genes_per_cluster_query are your lists of lists
# Calculate Overlap Index for each combination of clusters between the two lists of lists
overlap_results <- matrix(0, nrow = length(genes_per_cluster_ref), ncol = length(genes_per_cluster_query))
rownames(overlap_results) <- names(genes_per_cluster_ref)
colnames(overlap_results) <- names(genes_per_cluster_query)

for (i in seq_along(genes_per_cluster_ref)) {
  for (j in seq_along(genes_per_cluster_query)) {
    overlap_results[i, j] <- overlap_index(genes_per_cluster_ref[[i]], genes_per_cluster_query[[j]])
  }
}


# Convert the matrix to a long format suitable for ggplot
overlap_data_long <- as.data.frame(as.table(overlap_results))

# Renaming the columns for clarity
colnames(overlap_data_long) <- c("Cluster_ref", "Cluster_query", "OverlapIndex")

# Create the heatmap
heat_overlap_ref_query =    ggplot(overlap_data_long, aes(x = Cluster_ref, y = Cluster_query, fill = OverlapIndex)) +
  geom_tile(color="white") + # Use geom_tile() for heatmap squares
  scale_fill_viridis() + 
  theme_minimal() + # Use a minimal theme for cleaner appearance
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels for better readability
  labs(fill = "Overlap\nIndex", x = NULL, y = NULL, title = "Overlap Index Heatmap") # Set labels and title

print(heat_overlap_ref_query)
print(heat_jackard_ref_query + heat_overlap_ref_query )


pdf(file.path(PATH_ANALYSIS_OUTPUT_FIGURES, paste0("TOP",NUMBER_GENES_JACKARD,"Overlap.pdf")),  width = 6, height = 6)
print(heat_overlap_ref_query)
dev.off()

pdf(file.path(PATH_ANALYSIS_OUTPUT_FIGURES, paste0("TOP",NUMBER_GENES_JACKARD , "Jackard.pdf")),  width = 6, height = 6)
print(heat_jackard_ref_query)
dev.off()













