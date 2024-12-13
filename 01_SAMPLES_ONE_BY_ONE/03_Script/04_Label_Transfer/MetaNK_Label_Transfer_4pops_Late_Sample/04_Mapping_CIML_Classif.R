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


#Keep only relevant cells
if(DO_SUBSET_VARIABLE==TRUE){
  PBMC_Blood_All2 = SetIdent(PBMC_Blood_All , value=VARIABLE_TO_SUBSET )
  PBMC_Blood_All2 = subset(PBMC_Blood_All2, idents = CATEGORIES_TO_SUBSET)
}else{
  PBMC_Blood_All2 = PBMC_Blood_All
}



cat("\n\n")
print(p31)
cat("\n\n")

p32 = DimPlot(PBMC_Blood_All2, group.by = "predicted.id", pt.size = 0.6)

cat("\n\n")
print(p32)
cat("\n\n")

DimPlot(PBMC_Blood_All2 , group.by = "seurat_clusters")


cat("### Scoring on UMAP \n\n")

Markers_CIML = read.csv(file.path(PATH_EXPERIMENT_REFERENCE , "Tables_Gene_Signatures" , "Foltz_SI_2024", "Gene_Signatures.csv" ))

Markers_CIML %>%
  group_by(cluster) %>%
  #filter(!grepl("RPL|RPS|MT-", gene)) %>%
  top_n(n = FINDMARKERS_SHOWTOP, wt = avg_log2FC) -> top_CIML

Markers_CIML %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  #filter(!grepl("RPL|RPS|MT-", gene)) %>% #We use RPS and RPL for scoring but nor for plotind DotPlots
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_GENES_FOR_SCORING, wt = avg_log2FC) -> top_CIML


#Extracting the best genes 
list_per_cluster_ciml<- split(top_CIML, top_CIML$cluster)

# Extract the gene column from each data frame in the list
genes_per_cluster_ciml <- lapply(list_per_cluster_ciml, function(df) df$gene)

genes_per_cluster_ciml_restricted = lapply(genes_per_cluster_ciml, function(gene_list) {  intersect(gene_list, rownames(PBMC_Blood_All2))})

for (i in names(genes_per_cluster_ciml)){
  PBMC_Blood_All2 = AddModuleScore(PBMC_Blood_All2, features = list(genes_per_cluster_ciml_restricted[[i]]), pool= NULL ,name= i , seed=19)
}


#Ploting
features_to_plot <- paste0(names(genes_per_cluster_ciml), "1")
features_to_plot = gsub("-", ".", features_to_plot)
#features_to_plot = append(features_to_plot, "NKint1")
plot_titles <- c("CD56bright" ,"CD56dim" ,   "eML-1"   ,   "eML-2"   ,   "NKG2Cpos"  )

#plot_titles <- c("Score NK1", "Score NK3" , "Score NK2", "Score NKint")

# Generate the FeaturePlots with titles
plots <- lapply(seq_along(features_to_plot), function(i) {
  FeaturePlot(PBMC_Blood_All2, features = features_to_plot[i], pt.size = 0.6) &
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) &
    ggtitle(plot_titles[i])
})


# Combine the plots
combined_plot <- wrap_plots(plots)

# Display the combined plot
print(combined_plot)

pdf(file.path(PATH_ANALYSIS_OUTPUT_FIGURES, "Compare_Foltz" ,  "Scores_UMAP.pdf"),  width = 20, height = 4)
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
#All_Markers_ref = read.csv(file.path(PATH_EXPERIMENT_REFERENCE , "Tables_Gene_Signatures" , "Foltz_SI_2024", "Gene_Signatures.csv" ))

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
#PBMC_Blood_All = SetIdent(PBMC_Blood_All, value= "seurat_clusters")

#All_Markers_query = FindAllMarkers(PBMC_Blood_All2 , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

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

#Restrict to the genes that do appear in the query
genes_per_cluster_ref = lapply(genes_per_cluster_ref, function(gene_list) {  intersect(gene_list, rownames(PBMC_Blood_All2))})
genes_per_cluster_query = lapply(genes_per_cluster_query, function(gene_list) {  intersect(gene_list, rownames(PBMC_Blood_All2))})


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


pdf(file.path(PATH_ANALYSIS_OUTPUT_FIGURES, "Compare_Foltz", paste0("TOP",NUMBER_GENES_JACKARD,"Overlap.pdf")),  width = 6, height = 6)
print(heat_overlap_ref_query)
dev.off()

pdf(file.path(PATH_ANALYSIS_OUTPUT_FIGURES, "Compare_Foltz" , paste0("TOP",NUMBER_GENES_JACKARD , "Jackard.pdf")),  width = 6, height = 6)
print(heat_jackard_ref_query)
dev.off()













