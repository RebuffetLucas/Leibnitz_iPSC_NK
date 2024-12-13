## @knitr Prepare_Data

#Run a Full diagnostic for each sample
#Load the list of Seurat unfiltered objects

#Create the data list
LIST_RDS_UNFILTERED =  readRDS(file.path(PATH_EXPERIMENT_OUTPUT, "00_Create_List_Objects/List_Seurat_Objects_Before_QC.rds"))
LIST_NAMES_SAMPLES = names(LIST_RDS_UNFILTERED)

#Keep Only Metadata of interest

LIST_RDS_UNFILTERED_WITH_METADATA = LIST_RDS_UNFILTERED

# Add metadata to all individual dataframes
for (cell_id in LIST_NAMES_SAMPLES) {
  
  #Add pct ribo and pct mito
  LIST_RDS_UNFILTERED_WITH_METADATA[[cell_id]][["percent.mito"]] <- PercentageFeatureSet(LIST_RDS_UNFILTERED_WITH_METADATA[[cell_id]], pattern = "^MT-")
  LIST_RDS_UNFILTERED_WITH_METADATA[[cell_id]][['percent.ribo']] <- PercentageFeatureSet(object = LIST_RDS_UNFILTERED_WITH_METADATA[[cell_id]], pattern = "^RP[SL][[:digit:]]")
  
}



#Save the list of Seurat objects with Metadata if it does not already exist
if(file.exists(file.path(PATH_EXPERIMENT_OUTPUT, "00_Create_List_Objects/List_Seurat_Objects_Before_QC_With_Metadata.rds")) == FALSE){
  saveRDS(LIST_RDS_UNFILTERED_WITH_METADATA , file.path(PATH_EXPERIMENT_OUTPUT, "00_Create_List_Objects/List_Seurat_Objects_Before_QC_With_Metadata.rds"))
}


